import sys
import os
import re
import json
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool, Manager, Lock
from concurrent.futures import ProcessPoolExecutor, as_completed
from PyQt5.QtWidgets import (QApplication, QMainWindow, QFileDialog, QLabel, 
                            QVBoxLayout, QHBoxLayout, QWidget, QPushButton, 
                            QComboBox, QCheckBox, QGroupBox, QFormLayout,
                            QLineEdit, QMessageBox, QProgressBar, QTabWidget,
                            QDoubleSpinBox, QSpinBox, QRadioButton, QButtonGroup,
                            QScrollArea, QListWidget, QListWidgetItem, QDialog,
                            QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
                            QFrame, QSplitter, QInputDialog, QTextEdit, QSlider,QStatusBar,QSizePolicy)
from PyQt5.QtGui import QFont, QIcon, QColor, QIcon
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from collections import defaultdict
from Utils import try_bdf_file
from HelpDoc import get_topic_content

# Global functions for multiprocessing
def process_cbush_connection_chunk(args):
    """Process a chunk of CBUSH elements for connection identification"""
    cbush_chunk, model_data = args
    
    connections = {}
    processed_cbush = set()
    connection_id_start = cbush_chunk[0] if cbush_chunk else 1
    
    for i, cbush_id in enumerate(cbush_chunk):
        if cbush_id in processed_cbush:
            continue
            
        # Find ALL connected CBUSHes (including potential Airbus)
        all_connected_cbushes = find_all_connected_cbushes(cbush_id, model_data)
        
        if all_connected_cbushes:
            connection_id = connection_id_start + i
            connection = Connection(connection_id)
            
            # Add all CBUSHes to connection first
            for connected_cbush_id in all_connected_cbushes:
                if connected_cbush_id not in processed_cbush:
                    plates = get_cbush_plates(connected_cbush_id, model_data)
                    if len(plates) >= 2:
                        plate1_id, plate2_id = plates[0], plates[1]
                        connection.add_cbush(connected_cbush_id, plate1_id, plate2_id)
                        processed_cbush.add(connected_cbush_id)
            
            # Finalize connection with Airbus detection
            if connection.cbush_elements:
                connection.finalize_with_airbus_detection(model_data)
                connections[connection_id] = connection
    
    return connections, processed_cbush

def find_all_connected_cbushes(start_cbush_id, model_data):
    """Find ALL CBUSHes connected by shared nodes (including Airbus)"""
    connection_cbushes = []
    processed = set()
    stack = [start_cbush_id]
    
    while stack:
        current_cbush_id = stack.pop()
        if current_cbush_id in processed:
            continue
            
        connection_cbushes.append(current_cbush_id)
        processed.add(current_cbush_id)
        
        current_elem = model_data['elements'][current_cbush_id]
        
        # Find ALL other CBUSHes sharing nodes
        for node_id in current_elem.nodes:
            for other_elem_id, other_elem in model_data['elements'].items():
                if (other_elem.type == 'CBUSH' and 
                    other_elem_id not in processed and
                    node_id in other_elem.nodes):
                    stack.append(other_elem_id)
    
    return connection_cbushes

def get_cbush_plates(cbush_id, model_data):
    """Get the two plates connected by a CBUSH element"""
    cbush_elem = model_data['elements'][cbush_id]
    node1, node2 = cbush_elem.nodes
    
    # Get connected plates for both nodes
    plates_node1 = find_connected_plates(node1, model_data)
    plates_node2 = find_connected_plates(node2, model_data)
    
    # Return the first plate from each end
    plates = []
    if plates_node1:
        plates.append(plates_node1[0])
    if plates_node2:
        plates.append(plates_node2[0])
    
    return plates

def find_connected_plates(node_id, model_data):
    """Find plates connected to a node through RBE3"""
    connected_plates = []
    
    # Check if this node is dependent in any RBE3
    rbe3_independent_nodes = []
    for elem_id, elem in model_data['rigid_elements'].items():
        if hasattr(elem, 'refgrid') and elem.refgrid == node_id:
            if hasattr(elem, "independent_nodes") and len(elem.independent_nodes) > 0:
                if elem.independent_nodes[0]==node_id and len(elem.independent_nodes) > 1:
                    rbe3_independent_nodes.append(elem.independent_nodes[1])
                else:
                    rbe3_independent_nodes.append(elem.independent_nodes[0])
            break
    
    # Find plates connected to independent nodes
    if rbe3_independent_nodes:
        plates = find_direct_connected_plates(rbe3_independent_nodes[0], model_data)
        connected_plates.extend(plates)
    else:
        plates = find_direct_connected_plates(node_id, model_data)
        connected_plates.extend(plates)
    
    return connected_plates

def find_direct_connected_plates(node_id, model_data):
    """Find plates directly connected to a node"""
    connected_plates = []
    
    for elem_id, elem in model_data['elements'].items():
        if elem.type in ['CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6']:
            if node_id in elem.node_ids:
                connected_plates.append(elem_id)
                break
    
    return connected_plates

def calculate_stiffness_chunk(args):
    """Calculate stiffness for a chunk of elements"""
    element_chunk, model_data, connections_data, parameters_data, specs_dict = args
    
    calculator = HuthCalculator()
    stiffness_results = {}
    
    for elem_id, connection_id, params in element_chunk:
        # Set calculator parameters
        calculator.set_parameters({
            'fastener_e': params.get('fastener_e', 110.0),
            'fastener_type': params.get('fastener_type', 'bolt')
        })
        
        # Get connection
        connection = connections_data.get(connection_id)
        
        # Calculate stiffness
        stiffness = calculator.calculate_huth_stiffness(
            connection, elem_id, 
            params.get('fastener_diameter_mm', 6.35),
            params.get('connection_method', 'normal'),
            params.get('k4', 10.0),
            params.get('k5', 1e8),
            params.get('k6', 1e8),
            params.get('use_airbus_method', False),
            params.get('use_material_orientation', False)
        )
        
        stiffness_results[elem_id] = stiffness
    
    return stiffness_results

class Connection:
    """Represents a multi-plate connection with grouped CBUSH elements"""
    
    def __init__(self, connection_id):
        self.connection_id = connection_id
        self.cbush_elements = []  # ALL CBUSH element IDs in this connection
        self.regular_cbushes = []  # Only regular (non-Airbus) CBUSHes
        self.airbus_cbushes = []   # Airbus CBUSHes
        self.plates = []  # List of plate element IDs
        self.cbush_to_plates = {}  # Maps CBUSH ID to [plate1_id, plate2_id]
        self.plate_properties = {}  # Maps plate_id to properties dict
        self.total_thickness = 0.0
        self.ordered_plates = []  # Plates ordered from one end to the other
        self.airbus_cbush = None  # Primary Airbus CBUSH element ID
        self.is_airbus_method = False  # Flag to indicate if this connection uses Airbus method
        self.end_nodes = []  # End nodes for Airbus method
    
    def add_cbush(self, cbush_id, plate1_id, plate2_id):
        """Add a CBUSH element to this connection"""
        self.cbush_elements.append(cbush_id)
        self.cbush_to_plates[cbush_id] = [plate1_id, plate2_id]
        
        # Add plates to the connection if not already present
        if plate1_id not in self.plates:
            self.plates.append(plate1_id)
        if plate2_id not in self.plates:
            self.plates.append(plate2_id)
    
    def finalize_with_airbus_detection(self, model_data):
        """Finalize connection with proper Airbus detection"""
        # First, classify CBUSHes as regular vs Airbus
        self._classify_cbushes(model_data)
        
        # Order plates using only regular CBUSHes
        self._order_plates_by_regular_cbushes(model_data)
        
        # Calculate total thickness
        self.total_thickness = sum(
            self.plate_properties.get(plate_id, {}).get('thickness', 0.1) 
            for plate_id in self.ordered_plates
        )
        
        # Validate and set Airbus method
        if self.airbus_cbushes:
            self.is_airbus_method = True
            self.airbus_cbush = self.airbus_cbushes[0]  # Primary Airbus CBUSH
            self._find_end_nodes(model_data)
    
    def _classify_cbushes(self, model_data):
        """Classify CBUSHes as regular or Airbus based on stiffness properties"""
        for cbush_id in self.cbush_elements:
            if self._is_airbus_cbush(cbush_id, model_data):
                self.airbus_cbushes.append(cbush_id)
            else:
                self.regular_cbushes.append(cbush_id)
    
    def _is_airbus_cbush(self, cbush_id, model_data):
        """Check if a CBUSH has Airbus properties (high K1, low K2/K3)"""
        elem = model_data['elements'][cbush_id]
        if elem.pid in model_data['properties']:
            prop = model_data['properties'][elem.pid]
            if hasattr(prop, 'Ki') and len(prop.Ki) >= 3:
                k1, k2, k3 = prop.Ki[0], prop.Ki[1], prop.Ki[2]
                # Airbus method: high K1, low K2/K3
                return (k1 > 1000 and k2 < 1000 and k3 < 1000 and 
                       k1 > (max(k2, k3) * 10))
        return False
    
    def _order_plates_by_regular_cbushes(self, model_data):
        """Order plates using only regular CBUSHes"""
        if len(self.regular_cbushes) <= 1:
            self.ordered_plates = sorted(self.plates)
            return
        
        # Build node-to-regular-cbush mapping
        node_to_regular_cbush = defaultdict(list)
        for cbush_id in self.regular_cbushes:
            cbush_elem = model_data['elements'][cbush_id]
            for node_id in cbush_elem.nodes:
                node_to_regular_cbush[node_id].append(cbush_id)
        
        # Find end nodes: connected to only 1 regular CBUSH
        end_nodes = []
        for node_id, regular_cbush_list in node_to_regular_cbush.items():
            if len(regular_cbush_list) == 1:
                end_nodes.append(node_id)
        
        if len(end_nodes) < 2:
            self.ordered_plates = sorted(self.plates)
            return
        
        # Build plate sequence starting from one end
        ordered_plates = []
        used_cbushes = set()
        
        # Start from first end node
        start_node = end_nodes[0]
        current_cbush = node_to_regular_cbush[start_node][0]
        current_node = start_node
        
        while current_cbush and current_cbush not in used_cbushes:
            used_cbushes.add(current_cbush)
            
            # Get plates for this CBUSH
            plate1_id, plate2_id = self.cbush_to_plates[current_cbush]
            
            if not ordered_plates:
                # First CBUSH - determine plate order based on start node
                plate_for_start = self._get_plate_for_node(current_node, [plate1_id, plate2_id], model_data)
                other_plate = plate2_id if plate_for_start == plate1_id else plate1_id
                ordered_plates = [plate_for_start, other_plate]
                
                # Move to the other node of this CBUSH
                cbush_elem = model_data['elements'][current_cbush]
                node1, node2 = cbush_elem.nodes
                current_node = node2 if current_node == node1 else node1
            else:
                # Add new plate
                if plate1_id not in ordered_plates:
                    ordered_plates.append(plate1_id)
                    new_plate = plate1_id
                else:
                    ordered_plates.append(plate2_id)
                    new_plate = plate2_id
                
                # Move to node connected to the new plate
                cbush_elem = model_data['elements'][current_cbush]
                node1, node2 = cbush_elem.nodes
                
                plate_for_node1 = self._get_plate_for_node(node1, [plate1_id, plate2_id], model_data)
                if plate_for_node1 == new_plate:
                    current_node = node1
                else:
                    current_node = node2
            
            # Find next regular CBUSH
            next_cbush = None
            for cbush_id in node_to_regular_cbush[current_node]:
                if cbush_id not in used_cbushes:
                    next_cbush = cbush_id
                    break
            
            current_cbush = next_cbush
        
        self.ordered_plates = ordered_plates if ordered_plates else sorted(self.plates)
    
    def _get_plate_for_node(self, node_id, candidate_plates, model_data):
        """Determine which plate a node belongs to through RBE3"""
        for elem_id, elem in model_data['rigid_elements'].items():
            if hasattr(elem, 'refgrid') and elem.refgrid == node_id:
                # This node is dependent node of RBE3
                if hasattr(elem, 'Gijs') and len(elem.Gijs) > 0:
                    if len(elem.Gijs[0]) > 0:
                        independent_node = elem.Gijs[0][0]
                        # Find which plate contains this independent node
                        for plate_id in candidate_plates:
                            if plate_id in model_data['elements']:
                                plate_elem = model_data['elements'][plate_id]
                                if hasattr(plate_elem, 'node_ids') and independent_node in plate_elem.node_ids:
                                    return plate_id
        return candidate_plates[0]  # Fallback
    
    def _find_end_nodes(self, model_data):
        """Find end nodes for Airbus method"""
        if len(self.ordered_plates) < 2:
            return
        
        end_nodes = []
        outer_plates = [self.ordered_plates[0], self.ordered_plates[-1]]
        
        # Find nodes from regular CBUSHes connected to outer plates
        for cbush_id in self.regular_cbushes:
            plate1_id, plate2_id = self.cbush_to_plates[cbush_id]
            if plate1_id in outer_plates or plate2_id in outer_plates:
                cbush_elem = model_data['elements'][cbush_id]
                for node_id in cbush_elem.nodes:
                    plate_for_node = self._get_plate_for_node(node_id, [plate1_id, plate2_id], model_data)
                    if plate_for_node in outer_plates:
                        end_nodes.append(node_id)
        
        self.end_nodes = list(set(end_nodes))[:2]  # Keep unique, max 2
    
    def get_outer_nodes_for_new_airbus(self, model_data):
        """Get nodes for creating new Airbus CBUSH"""
        if not self.end_nodes:
            self._find_end_nodes(model_data)
        return self.end_nodes[:2]

class HuthCalculator:
    """Class for calculating Huth stiffness for CBUSH elements"""
    
    def __init__(self):
        self.fastener_e = 110.0
        self.fastener_type = "bolt"
    
    def calculate_huth_stiffness(self, connection, cbush_id, fastener_diameter_mm, 
                                connection_method="normal", k4=10.0, k5=1e8, k6=1e8, 
                                use_airbus_method=False, use_material_orientation=False):
        """Calculate CBUSH stiffness"""
        fastener_diameter = fastener_diameter_mm
        fastener_e_mpa = self.fastener_e * 1000
        print(f"CBUSH {cbush_id}: Connection ID {connection.connection_id if connection else 'None'}")
        print(f"CBUSH {cbush_id}: In cbush_to_plates: {cbush_id in connection.cbush_to_plates if connection else 'N/A'}")
        if not connection:
            return {
                'K1': 1e7, 'K2': 1e7, 'K3': 1e7,
                'K4': k4, 'K5': k5, 'K6': k6
            }
        
        # Check if this is an Airbus CBUSH
        is_airbus_cbush = (cbush_id in connection.airbus_cbushes)
        
        # For Airbus method
        if use_airbus_method and len(connection.plates) > 2:
            return self._calculate_airbus_stiffness(connection, cbush_id, fastener_diameter, 
                                                   fastener_e_mpa, k4, k5, k6, connection_method, 
                                                   use_material_orientation)
        
        # Handle existing Airbus CBUSHes even in non-Airbus mode
        if is_airbus_cbush:
            if connection.is_airbus_method:
                # Keep existing Airbus stiffness distribution
                total_thickness = connection.total_thickness
                k_axial = self._calculate_axial_stiffness_total(total_thickness, fastener_diameter, fastener_e_mpa)
                return {
                    'K1': k_axial,
                    'K2': 100.0,
                    'K3': 100.0,
                    'K4': k4, 'K5': k5, 'K6': k6
                }
            else:
                # Convert existing Airbus CBUSH to normal
                if cbush_id in connection.cbush_to_plates:
                    cbush_plates = connection.cbush_to_plates[cbush_id]
                    plate1_id, plate2_id = cbush_plates
                    plate1_props = connection.plate_properties.get(plate1_id, self._get_default_plate_props())
                    plate2_props = connection.plate_properties.get(plate2_id, self._get_default_plate_props())
                    
                    shear_result = self._calculate_huth_stiffness_document(
                        plate1_props, plate2_props, fastener_diameter, fastener_e_mpa, n=1, 
                        use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
                    
                    # Handle directional results
                    if isinstance(shear_result, tuple) and len(shear_result) == 2:
                        shear_stiffness_x, shear_stiffness_y = shear_result
                    else:
                        shear_stiffness_x = shear_stiffness_y = shear_result
                    
                    k_axial = self._calculate_axial_stiffness(plate1_props, plate2_props, 
                                                              fastener_diameter, fastener_e_mpa)
                    
                    return {
                        'K1': k_axial,
                        'K2': shear_stiffness_x,
                        'K3': shear_stiffness_y,
                        'K4': k4, 'K5': k5, 'K6': k6
                    }
        
        # Regular CBUSH calculation
        if cbush_id not in connection.cbush_to_plates:
            return {
                'K1': 1e7, 'K2': 1e7, 'K3': 1e7,
                'K4': k4, 'K5': k5, 'K6': k6
            }
        
        cbush_plates = connection.cbush_to_plates[cbush_id]
        plate1_id, plate2_id = cbush_plates
        
        plate1_props = connection.plate_properties.get(plate1_id, self._get_default_plate_props())
        plate2_props = connection.plate_properties.get(plate2_id, self._get_default_plate_props())
        
        # Calculate stiffness based on method
        if connection_method == "two_plate":
            shear_result = self._calculate_huth_stiffness_document(
                plate1_props, plate2_props, fastener_diameter, fastener_e_mpa, n=1,
                use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
        elif connection_method == "hypermesh":
            shear_result = self._calculate_hypermesh_huth_stiffness(
                connection, cbush_id, fastener_diameter, fastener_e_mpa, use_material_orientation)
        else:
            shear_result = self._calculate_multiplate_huth_stiffness(
                connection, cbush_id, fastener_diameter, fastener_e_mpa, use_material_orientation)
        
        # Handle directional results
        if isinstance(shear_result, tuple) and len(shear_result) == 2:
            shear_stiffness_x, shear_stiffness_y = shear_result
        else:
            shear_stiffness_x = shear_stiffness_y = shear_result
        
        k_axial = self._calculate_axial_stiffness(plate1_props, plate2_props, 
                                                  fastener_diameter, fastener_e_mpa)
        
        return {
            'K1': k_axial,
            'K2': shear_stiffness_x,
            'K3': shear_stiffness_y,
            'K4': k4, 'K5': k5, 'K6': k6
        }
    
    def _calculate_airbus_stiffness(self, connection, cbush_id, d, E_f, k4, k5, k6, 
                                   connection_method, use_material_orientation):
        """Calculate stiffness for Airbus method with connection method awareness"""
        is_airbus_cbush = (cbush_id in connection.airbus_cbushes)
        
        if is_airbus_cbush:
            # Airbus CBUSH - only axial stiffness
            total_thickness = connection.total_thickness
            k_axial = self._calculate_axial_stiffness_total(total_thickness, d, E_f)
            
            return {
                'K1': k_axial,
                'K2': 100.0,
                'K3': 100.0,
                'K4': k4, 'K5': k5, 'K6': k6
            }
        else:
            # Regular CBUSH in Airbus method - calculate shear stiffness using selected method
            if cbush_id in connection.cbush_to_plates:
                cbush_plates = connection.cbush_to_plates[cbush_id]
                plate1_id, plate2_id = cbush_plates
                plate1_props = connection.plate_properties.get(plate1_id, self._get_default_plate_props())
                plate2_props = connection.plate_properties.get(plate2_id, self._get_default_plate_props())
                
                # Use the selected connection method for shear stiffness calculation
                if connection_method == "two_plate":
                    shear_result = self._calculate_huth_stiffness_document(
                        plate1_props, plate2_props, d, E_f, n=1,
                        use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
                elif connection_method == "hypermesh":
                    shear_result = self._calculate_hypermesh_huth_stiffness(
                        connection, cbush_id, d, E_f, use_material_orientation)
                else:  # normal
                    shear_result = self._calculate_multiplate_huth_stiffness(
                        connection, cbush_id, d, E_f, use_material_orientation)
            else:
                shear_result = 1e7
            
            # Handle directional results
            if isinstance(shear_result, tuple) and len(shear_result) == 2:
                shear_stiffness_x, shear_stiffness_y = shear_result
            else:
                shear_stiffness_x = shear_stiffness_y = shear_result
            
            return {
                'K1': 100.0,
                'K2': shear_stiffness_x,
                'K3': shear_stiffness_y,
                'K4': k4, 'K5': k5, 'K6': k6
            }
    
    def _calculate_axial_stiffness_total(self, total_thickness, d, E_f):
        """Calculate axial stiffness using total connection thickness"""
        pi = np.pi
        k_axial = (E_f * pi * d**2) / (4 * total_thickness) if total_thickness > 0 else 1e12
        return k_axial
    
    def _calculate_hypermesh_huth_stiffness(self, connection, cbush_id, d, E_f, use_material_orientation):
        """Calculate HyperMesh-style multi-plate Huth stiffness"""
        num_plates = len(connection.plates)
        
        if num_plates <= 2:
            cbush_plates = connection.cbush_to_plates[cbush_id]
            plate1_props = connection.plate_properties.get(cbush_plates[0], self._get_default_plate_props())
            plate2_props = connection.plate_properties.get(cbush_plates[1], self._get_default_plate_props())
            return self._calculate_huth_stiffness_document(
                plate1_props, plate2_props, d, E_f, n=1,
                use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
        
        # Multi-plate connection
        cbush_plates = connection.cbush_to_plates[cbush_id]
        plate1_id, plate2_id = cbush_plates
        ordered_plates = connection.ordered_plates
        
        plate1_is_outer = (plate1_id == ordered_plates[0] or plate1_id == ordered_plates[-1])
        plate2_is_outer = (plate2_id == ordered_plates[0] or plate2_id == ordered_plates[-1])
        
        plate1_props = connection.plate_properties.get(plate1_id, self._get_default_plate_props())
        plate2_props = connection.plate_properties.get(plate2_id, self._get_default_plate_props())
        
        if plate1_is_outer or plate2_is_outer:
            if plate1_is_outer:
                outer_props = plate1_props
                middle_props = plate2_props
            else:
                outer_props = plate2_props
                middle_props = plate1_props
        else:
            strength1 = self._calculate_plate_strength(plate1_props)
            strength2 = self._calculate_plate_strength(plate2_props)
            
            if strength1 >= strength2:
                middle_props = plate1_props
                outer_props = plate2_props
            else:
                middle_props = plate2_props
                outer_props = plate1_props
        
        return self._calculate_huth_stiffness_document(
            middle_props, outer_props, d, E_f, n=2,
            use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
    
    def _calculate_multiplate_huth_stiffness(self, connection, cbush_id, d, E_f, use_material_orientation):
        """Calculate multi-plate Huth stiffness"""
        num_plates = len(connection.plates)
        
        if num_plates <= 2:
            cbush_plates = connection.cbush_to_plates[cbush_id]
            plate1_props = connection.plate_properties.get(cbush_plates[0], self._get_default_plate_props())
            plate2_props = connection.plate_properties.get(cbush_plates[1], self._get_default_plate_props())
            return self._calculate_huth_stiffness_document(
                plate1_props, plate2_props, d, E_f, n=1,
                use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
        
        elif num_plates == 3:
            ordered_plates = connection.ordered_plates
            if len(ordered_plates) >= 3:
                outer_plate1_props = connection.plate_properties.get(ordered_plates[0], self._get_default_plate_props())
                middle_plate_props = connection.plate_properties.get(ordered_plates[1], self._get_default_plate_props())
                outer_plate2_props = connection.plate_properties.get(ordered_plates[2], self._get_default_plate_props())
                
                strength1 = self._calculate_plate_strength(outer_plate1_props)
                strength2 = self._calculate_plate_strength(outer_plate2_props)
                
                if strength1 <= strength2:
                    outer_equiv_props = outer_plate1_props.copy()
                else:
                    outer_equiv_props = outer_plate2_props.copy()
                
                return self._calculate_huth_stiffness_document(
                    middle_plate_props, outer_equiv_props, d, E_f, n=2,
                    use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
            else:
                cbush_plates = connection.cbush_to_plates[cbush_id]
                plate1_props = connection.plate_properties.get(cbush_plates[0], self._get_default_plate_props())
                plate2_props = connection.plate_properties.get(cbush_plates[1], self._get_default_plate_props())
                return self._calculate_huth_stiffness_document(
                    plate1_props, plate2_props, d, E_f, n=1,
                    use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
        
        else:
            ordered_plates = connection.ordered_plates
            plate_props_list = [connection.plate_properties.get(pid, self._get_default_plate_props()) 
                              for pid in ordered_plates]
            
            strongest_idx = 0
            max_strength = 0
            for i, props in enumerate(plate_props_list):
                strength = self._calculate_plate_strength(props)
                if strength > max_strength:
                    max_strength = strength
                    strongest_idx = i
            
            middle_plate_props = plate_props_list[strongest_idx]
            
            other_plates = [props for i, props in enumerate(plate_props_list) if i != strongest_idx]
            if other_plates:
                avg_thickness = sum(props.get('thickness', 1.0) for props in other_plates) / len(other_plates)
                avg_modulus = sum(props.get('elastic_modulus', 70000.0) for props in other_plates) / len(other_plates)
                
                outer_equiv_props = {
                    'type': other_plates[0].get('type', 'metal'),
                    'thickness': avg_thickness,
                    'elastic_modulus': avg_modulus
                }
            else:
                outer_equiv_props = self._get_default_plate_props()
            
            return self._calculate_huth_stiffness_document(
                middle_plate_props, outer_equiv_props, d, E_f, n=2,
                use_material_orientation=use_material_orientation, cbush_id=cbush_id, connection=connection)
    
    def _calculate_plate_strength(self, plate_props):
        """Calculate relative strength measure"""
        thickness = plate_props.get('thickness', 1.0)
        modulus = plate_props.get('elastic_modulus', 70000.0)
        strength = thickness * np.sqrt(modulus / 1000.0)
        return strength
    
    def _get_default_plate_props(self):
        """Get default plate properties"""
        return {
            'type': 'metal',
            'thickness': 1.0,
            'elastic_modulus': 70000.0,
            'elastic_modulus_x': 70000.0,
            'elastic_modulus_y': 70000.0,
            'use_directional': False
        }
    
    def _calculate_huth_stiffness_document(self, plate1_props, plate2_props, d, E_f, n, 
                                          use_material_orientation=False, cbush_id=None, connection=None):
        """Calculate Huth stiffness according to document formulation with material orientation"""
        
        if not use_material_orientation or cbush_id is None or connection is None:
            # Original effective modulus approach
            t1 = plate1_props.get('thickness', 1.0)
            t2 = plate2_props.get('thickness', 1.0)
            E1 = plate1_props.get('elastic_modulus', 70000.0)
            E2 = plate2_props.get('elastic_modulus', 70000.0)

            a, b1, b2 = self._get_huth_parameters(plate1_props, plate2_props)
            
            b1 = b1 / n
            b2 = b2 / (n ** 2)
            
            tmp0 = (t1 + t2) / (2 * d)
            tmp1 = 1/(t1*E1) + 1/(2*t1*E_f) if t1 > 0 else 0
            tmp2 = 1/(t2*E2) + 1/(2*t2*E_f) if t2 > 0 else 0
            
            compliance = (tmp0**a) * b1 * tmp1 + (tmp0**a) * b2 * tmp2
            stiffness = 1.0 / compliance if compliance > 0 else 1e12
            
            return stiffness
        
        else:
            # Material orientation approach - calculate directional stiffnesses
            return self._calculate_directional_huth_stiffness(
                connection, cbush_id, plate1_props, plate2_props, d, E_f, n)
    
    def _calculate_directional_huth_stiffness(self, connection, cbush_id, plate1_props, plate2_props, d, E_f, n):
        """Unified directional Huth stiffness calculation"""
        
        model_data = getattr(connection, '_model_data', None)
        if not model_data or cbush_id not in model_data.get('elements', {}):
            return self._calculate_huth_stiffness_document(
                plate1_props, plate2_props, d, E_f, n, use_material_orientation=False)
        
        cbush_elem = model_data['elements'][cbush_id]
        
        # Step 1: Get fastener coordinate system (handles both aligned and non-aligned)
        fastener_system = self._get_fastener_coordinate_system(cbush_elem, model_data)
        
        # Step 2: Get plate elements
        plate1_id, plate2_id = connection.cbush_to_plates[cbush_id]
        plate1_elem = model_data['elements'].get(plate1_id)
        plate2_elem = model_data['elements'].get(plate2_id)
        
        if not plate1_elem or not plate2_elem:
            return self._calculate_huth_stiffness_document(
                plate1_props, plate2_props, d, E_f, n, use_material_orientation=False)
        
        # Step 3: Calculate effective angles for each direction and each plate
        # This automatically handles both aligned and non-aligned cases
        theta1_e2, theta1_e3 = self._get_material_angles_relative_to_fastener(
            fastener_system, plate1_elem, plate1_props, model_data)
        theta2_e2, theta2_e3 = self._get_material_angles_relative_to_fastener(
            fastener_system, plate2_elem, plate2_props, model_data)
        
        # Step 4: Transform material properties for each plate and direction
        E1_e2, E1_e3 = self._get_directional_moduli(plate1_props, theta1_e2, theta1_e3)
        E2_e2, E2_e3 = self._get_directional_moduli(plate2_props, theta2_e2, theta2_e3)
        print(f"theta1_e2: {np.degrees(theta1_e2):.1f}°")
        print(f"theta1_e3: {np.degrees(theta1_e3):.1f}°")
        print(f"theta2_e2: {np.degrees(theta2_e2):.1f}°")
        print(f"theta2_e3: {np.degrees(theta2_e3):.1f}°")
        print(f"Transformed E1_e2: {E1_e2:.1f}, E1_e3: {E1_e3:.1f}")
        print(f"Transformed E2_e2: {E2_e2:.1f}, E2_e3: {E2_e3:.1f}")
        # Step 5: Calculate Huth stiffness
        return self._calculate_huth_with_directional_moduli(
            E1_e2, E1_e3, E2_e2, E2_e3, plate1_props, plate2_props, d, E_f, n)

    def _calculate_huth_with_directional_moduli(self, E1_e2, E1_e3, E2_e2, E2_e3, 
                                            plate1_props, plate2_props, d, E_f, n):
        """Calculate Huth stiffness using directional moduli"""
        
        t1 = plate1_props.get('thickness', 1.0)
        t2 = plate2_props.get('thickness', 1.0)
        
        a, b1, b2 = self._get_huth_parameters(plate1_props, plate2_props)
        b1 = b1 / n
        b2 = b2 / (n ** 2)
        
        tmp0 = (t1 + t2) / (2 * d)
        
        # Calculate K2 (fastener e2 direction)
        tmp1_e2 = 1/(t1*E1_e2) + 1/(2*t1*E_f) if t1 > 0 and E1_e2 > 0 else 0
        tmp2_e2 = 1/(t2*E2_e2) + 1/(2*t2*E_f) if t2 > 0 and E2_e2 > 0 else 0
        compliance_e2 = (tmp0**a) * b1 * tmp1_e2 + (tmp0**a) * b2 * tmp2_e2
        K2 = 1.0 / compliance_e2 if compliance_e2 > 0 else 1e12
        
        # Calculate K3 (fastener e3 direction)
        tmp1_e3 = 1/(t1*E1_e3) + 1/(2*t1*E_f) if t1 > 0 and E1_e3 > 0 else 0
        tmp2_e3 = 1/(t2*E2_e3) + 1/(2*t2*E_f) if t2 > 0 and E2_e3 > 0 else 0
        compliance_e3 = (tmp0**a) * b1 * tmp1_e3 + (tmp0**a) * b2 * tmp2_e3
        K3 = 1.0 / compliance_e3 if compliance_e3 > 0 else 1e12
        
        return K2, K3

    def _get_material_angles_relative_to_fastener(self, fastener_system, plate_elem, plate_props, model_data):
        """Calculate angle between material direction and fastener directions"""
        
        # Step 1: Get plate geometry
        plate_normal = self._get_plate_normal(plate_elem, model_data)
        elem_x_direction = self._get_corrected_element_x_direction(plate_elem, model_data)
        
        # Step 2: Get material direction in global coordinates  
        material_dir_global = self._get_material_direction_global(
            plate_elem, elem_x_direction, plate_normal, model_data)
        
        # Step 3: Project fastener directions onto plate surface
        fastener_e2_on_plate = self._project_vector_onto_plane(
            fastener_system['e2'], plate_normal)
        fastener_e3_on_plate = self._project_vector_onto_plane(
            fastener_system['e3'], plate_normal)
        
        # Step 4: Calculate angles 
        theta_e2 = self._angle_between_vectors(material_dir_global, fastener_e2_on_plate)
        theta_e3 = self._angle_between_vectors(material_dir_global, fastener_e3_on_plate)
        
        return theta_e2, theta_e3

    def _get_material_direction_global(self, plate_elem, elem_x_direction, plate_normal, model_data):
        """Get material direction in global coordinates for any theta_mcid"""
        
        try:
            theta_mcid = getattr(plate_elem, 'theta_mcid', 0)
            
            if theta_mcid == 0:
                # No rotation - material aligned with element x-direction
                return elem_x_direction
            
            elif isinstance(theta_mcid, (int, float)) and theta_mcid > 0 and isinstance(theta_mcid, int):
                # Coordinate system reference (MCID)
                return self._get_material_direction_from_coord_system(
                    theta_mcid, elem_x_direction, plate_normal, model_data)
            
            else:
                # Direct angle in degrees
                angle_rad = np.radians(float(theta_mcid))
                return self._rotate_vector_in_plane(elem_x_direction, plate_normal, angle_rad)
                
        except (AttributeError, TypeError, ValueError):
            return elem_x_direction

    def _get_material_direction_from_coord_system(self, mcid, elem_x_direction, plate_normal, model_data):
        """Handle coordinate system references (MCID)"""
        
        try:
            if mcid not in model_data.get('coords', {}):
                return elem_x_direction
            
            coord_sys = model_data['coords'][mcid]
            
            # Get coordinate system x-axis in global coordinates
            if hasattr(coord_sys, 'i'):
                cs_x_direction = np.array(coord_sys.i)
            elif hasattr(coord_sys, 'e1'):
                cs_x_direction = np.array(coord_sys.e1)
            else:
                return elem_x_direction
            
            # Project coordinate system x-direction onto plate surface
            cs_x_on_plate = self._project_vector_onto_plane(cs_x_direction, plate_normal)
            
            return cs_x_on_plate / np.linalg.norm(cs_x_on_plate) if np.linalg.norm(cs_x_on_plate) > 1e-10 else elem_x_direction
            
        except Exception:
            return elem_x_direction

    def _get_plate_normal(self, plate_elem, model_data):
        """Get plate normal vector correctly for any plate orientation"""
        
        try:
            # Get first 3 nodes of the plate
            if plate_elem.type in ['CQUAD4', 'CQUAD8']:
                node_ids = plate_elem.nodes[:3]  
            elif plate_elem.type in ['CTRIA3', 'CTRIA6']:
                node_ids = plate_elem.nodes[:3]
            else:
                return np.array([0.0, 0.0, 1.0])  
            
            # Get node positions
            positions = []
            for node_id in node_ids:
                if node_id in model_data.get('nodes', {}):
                    positions.append(np.array(model_data['nodes'][node_id].xyz))
                else:
                    return np.array([0.0, 0.0, 1.0])
            
            if len(positions) < 3:
                return np.array([0.0, 0.0, 1.0])
            
            # Calculate normal using cross product
            v1 = positions[1] - positions[0]  # G1 to G2
            v2 = positions[2] - positions[0]  # G1 to G3
            
            normal = np.cross(v1, v2)
            magnitude = np.linalg.norm(normal)
            
            if magnitude > 1e-10:
                return normal / magnitude
            else:
                return np.array([0.0, 0.0, 1.0])
                
        except Exception:
            return np.array([0.0, 0.0, 1.0])

    def _get_corrected_element_x_direction(self, plate_elem, model_data):
        """Get element x-direction correctly based on NASTRAN documentation"""
        
        try:
            if plate_elem.type in ['CTRIA3', 'CTRIA6']:
                # For triangular elements: x-axis is from G1 to G2 (side 1-2)
                n1_id, n2_id = plate_elem.nodes[0], plate_elem.nodes[1]
                
            elif plate_elem.type in ['CQUAD4', 'CQUAD8']:
                # For quad elements: coordinate system is more complex
                # According to NASTRAN docs, we need proper bisection calculation
                # For simplicity in fastener analysis, we'll use G1-G2 as approximation
                # but note this could be improved for high-precision applications
                n1_id, n2_id = plate_elem.nodes[0], plate_elem.nodes[1]
                
            else:
                return np.array([1.0, 0.0, 0.0])
            
            # Get node positions
            nodes = model_data.get('nodes', {})
            if n1_id in nodes and n2_id in nodes:
                n1_pos = np.array(nodes[n1_id].xyz)
                n2_pos = np.array(nodes[n2_id].xyz)
                
                x_direction = n2_pos - n1_pos
                magnitude = np.linalg.norm(x_direction)
                
                if magnitude > 1e-10:
                    return x_direction / magnitude
                else:
                    return np.array([1.0, 0.0, 0.0])
            else:
                return np.array([1.0, 0.0, 0.0])
                
        except Exception:
            return np.array([1.0, 0.0, 0.0])

    def _project_vector_onto_plane(self, vector, normal):
        """Project any vector onto any plane defined by normal"""
        
        try:
            vector = np.array(vector)
            normal = np.array(normal)
            
            # Normalize normal
            normal_magnitude = np.linalg.norm(normal)
            if normal_magnitude < 1e-10:
                return vector
            normal = normal / normal_magnitude
            
            # Project: v_proj = v - (v·n)n
            dot_product = np.dot(vector, normal)
            projection_on_normal = dot_product * normal
            projected_vector = vector - projection_on_normal
            
            # Normalize
            magnitude = np.linalg.norm(projected_vector)
            if magnitude > 1e-10:
                return projected_vector / magnitude
            else:
                # If vector is parallel to normal, return arbitrary perpendicular vector
                return self._get_perpendicular_vector(normal)
                
        except Exception:
            return np.array([1.0, 0.0, 0.0])

    def _rotate_vector_in_plane(self, vector, normal, angle_rad):
        """Rotate vector in plane using Rodrigues' rotation formula"""
        
        try:
            vector = np.array(vector)
            normal = np.array(normal)
            
            # Normalize
            vector = vector / np.linalg.norm(vector)
            normal = normal / np.linalg.norm(normal)
            
            # Rodrigues' rotation formula: v' = v*cos(θ) + (n×v)*sin(θ) + n*(n·v)*(1-cos(θ))
            cos_angle = np.cos(angle_rad)
            sin_angle = np.sin(angle_rad)
            
            cross_product = np.cross(normal, vector)
            dot_product = np.dot(normal, vector)
            
            rotated = (vector * cos_angle + 
                    cross_product * sin_angle + 
                    normal * dot_product * (1 - cos_angle))
            
            return rotated / np.linalg.norm(rotated)
            
        except Exception:
            return vector

    def _angle_between_vectors(self, v1, v2):
        """Calculate angle between any two vectors"""
        
        try:
            v1 = np.array(v1)
            v2 = np.array(v2)
            
            # Normalize
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            
            # Calculate angle
            cos_angle = np.dot(v1, v2)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            
            return np.arccos(abs(cos_angle))  # Always positive angle
            
        except Exception:
            return 0.0

    def _get_perpendicular_vector(self, vector):
        """Get any vector perpendicular to given vector"""
        
        vector = np.array(vector)
        
        # Find component with smallest absolute value
        abs_components = np.abs(vector)
        min_idx = np.argmin(abs_components)
        
        # Create perpendicular vector
        perp = np.zeros(3)
        perp[min_idx] = 1.0
        
        # Make it perpendicular: perp = perp - (perp·v)v
        dot_product = np.dot(perp, vector)
        perp = perp - dot_product * vector
        
        return perp / np.linalg.norm(perp)

    def _get_directional_moduli(self, plate_props, theta_e2, theta_e3):
        """Get moduli with proper isotropic handling"""
        
        # Check if isotropic first
        plate_type = plate_props.get('type', 'metal')
        if plate_type == 'metal' or 'isotropic' in str(plate_type).lower():
            E_iso = plate_props.get('elastic_modulus', 70000.0)
            return E_iso, E_iso  # Same for both directions
        
        # Orthotropic/composite case
        if 'E1' in plate_props and 'E2' in plate_props:
            E1, E2 = plate_props['E1'], plate_props['E2']
            G12 = plate_props.get('G12', 26900.0)
            nu12 = plate_props.get('NU12', 0.3)
        else:
            E_eff = plate_props.get('elastic_modulus', 70000.0)
            return E_eff, E_eff
        
        # Transform for each direction using simple Jones formula
        E_e2 = self._jones_simple_transform(E1, E2, G12, nu12, theta_e2)
        E_e3 = self._jones_simple_transform(E1, E2, G12, nu12, theta_e3)
        
        return E_e2, E_e3
    
    def _get_fastener_coordinate_system(self, cbush_elem, model_data):
        """Get fastener coordinate system like TCL fastenerSystem function"""
        
        # Get fastener nodes
        node1_id, node2_id = cbush_elem.nodes[:2]
        node1_pos = np.array(model_data['nodes'][node1_id].xyz) if node1_id in model_data['nodes'] else np.array([0,0,0])
        node2_pos = np.array(model_data['nodes'][node2_id].xyz) if node2_id in model_data['nodes'] else np.array([1,0,0])
        
        # Fastener e1 axis (along fastener)
        e1 = node2_pos - node1_pos
        e1_magnitude = np.linalg.norm(e1)
        if e1_magnitude > 1e-10:
            e1 = e1 / e1_magnitude
        else:
            e1 = np.array([1.0, 0.0, 0.0])
        
        # Get orientation vector from CBUSH
        vector = self._get_cbush_orientation_vector(cbush_elem)
        
        # Calculate e3 (like TCL code)
        if np.linalg.norm(vector) > 1e-10:
            e3 = np.cross(e1, vector)
            e3_magnitude = np.linalg.norm(e3)
            if e3_magnitude > 1e-10:
                e3 = e3 / e3_magnitude
                e2 = np.cross(e3, e1)
            else:
                e2, e3 = self._get_default_perpendicular_axes(e1)
        else:
            e2, e3 = self._get_default_perpendicular_axes(e1)
        
        return {
            'e1': e1,  # Along fastener
            'e2': e2,  # Fastener y-direction (used for K2)
            'e3': e3   # Fastener z-direction (used for K3)
        }
    
    def _get_cbush_orientation_vector(self, cbush_elem):
        """Get orientation vector from CBUSH element"""
        try:
            # Try to get orientation vector from CBUSH element
            if hasattr(cbush_elem, 'x') and cbush_elem.x is not None:
                if isinstance(cbush_elem.x, (list, tuple, np.ndarray)) and len(cbush_elem.x) >= 3:
                    return np.array(cbush_elem.x[:3])
                    
            # Try to get from g0 (orientation node)
            if hasattr(cbush_elem, 'g0') and cbush_elem.g0 is not None and cbush_elem.g0 > 0:
                # This would require access to model_data['nodes']
                return np.array([0.0, 1.0, 0.0])  # Default for now
                    
            # Default vector
            return np.array([0.0, 0.0, 0.0])
            
        except (AttributeError, IndexError, KeyError):
            return np.array([0.0, 0.0, 0.0])
    
    def _get_default_perpendicular_axes(self, e1):
        """Get perpendicular axes when no orientation vector available (like TCL MCID=-1 case)"""
        # Find the component with smallest absolute value
        abs_components = np.abs(e1)
        min_idx = np.argmin(abs_components)
        
        # Create unit vector in that direction
        b = np.zeros(3)
        b[min_idx] = 1.0
        
        # Generate e2 perpendicular to e1
        dot_product = np.dot(e1, b)
        projection = dot_product * e1
        e2 = b - projection
        e2 = e2 / np.linalg.norm(e2)
        
        # Generate e3 = e1 × e2
        e3 = np.cross(e1, e2)
        
        return e2, e3
    
    def _jones_simple_transform(self, E1, E2, G12, nu12, theta_rad):
        """Apply simple Jones transformation formula (directional, not extended)"""
        try:
            cos_th = np.cos(theta_rad)
            sin_th = np.sin(theta_rad)
            
            # Simple directional Jones formula 
            inv_E = ((1.0/E1) * cos_th**4 + 
                    (1.0/E2) * sin_th**4 + 
                    (1.0/G12 - 2.0*nu12/E1) * cos_th**2 * sin_th**2)
            
            E_transformed = 1.0 / inv_E if inv_E > 0 else E1
            return E_transformed
            
        except (ZeroDivisionError, ValueError):
            return (E1 * E2)**0.5
    
    def _get_huth_parameters(self, plate1_props, plate2_props):
        """Get Huth parameters based on fastener and plate types"""
        a = 2/3
        b1 = 3.0
        b2 = 3.0
        
        if "rivet" in self.fastener_type.lower():
            a = 2/5
            b1 = 2.2
            b2 = 2.2
        
        if "composite" in plate1_props.get('type', 'metal').lower():
            b1 = 4.2
        if "composite" in plate2_props.get('type', 'metal').lower():
            b2 = 4.2
        
        return a, b1, b2
    
    def _calculate_axial_stiffness(self, plate1_props, plate2_props, d, E_f):
        """Calculate axial stiffness for fastener"""
        t1 = plate1_props.get('thickness', 1.0)
        t2 = plate2_props.get('thickness', 1.0)
        total_thickness = t1 + t2
        
        pi = np.pi
        k_axial = (E_f * pi * d**2) / (4 * total_thickness) if total_thickness > 0 else 1e12
        
        return k_axial
    
    def set_parameters(self, parameters):
        """Set parameters"""
        self.fastener_e = parameters.get('fastener_e', 110.0)
        self.fastener_type = parameters.get('fastener_type', 'bolt')

class FastenerSpec:
    """Class to store fastener specification details"""
    
    def __init__(self, spec_name=""):
        self.spec_name = spec_name
        self.fastener_type = "Bolt"
        self.material = "Titanium"
        self.elastic_modulus = 110.0
    
    def to_dict(self):
        return {
            'spec_name': self.spec_name,
            'fastener_type': self.fastener_type,
            'material': self.material,
            'elastic_modulus': self.elastic_modulus
        }
    
    @classmethod
    def from_dict(cls, data):
        spec = cls(data.get('spec_name', ''))
        spec.fastener_type = data.get('fastener_type', 'Bolt')
        spec.material = data.get('material', 'Titanium')
        spec.elastic_modulus = data.get('elastic_modulus', 110.0)
        return spec

class FastenerGroup:
    """Class to store group of elements with specific fastener configuration"""
    
    def __init__(self, set_id=0, set_name="", spec_name="", diameter_mm=6.35):
        self.set_id = set_id
        self.set_name = set_name
        self.original_name = set_name
        self.spec_name = spec_name
        self.diameter_mm = diameter_mm
        self.element_ids = []
        self.manual = False
        self.problematic = False
        self.diameter_missing = False
        self.export_name = set_name
    
    def update_set_name(self):
        if self.manual:
            return
        
        parts = self.set_name.split('__')
        if len(parts) >= 4:
            parts[2] = f"{self.diameter_mm:.2f}"
            if self.spec_name:
                parts[3] = self.spec_name
            self.set_name = '__'.join(parts)
    
    def update_export_name(self):
        if not self.original_name:
            self.export_name = self.set_name
            return
            
        parts = self.original_name.split('__')
        if len(parts) >= 4:
            if self.spec_name:
                parts[3] = self.spec_name
            parts[2] = f"{self.diameter_mm:.2f}"
            self.export_name = '__'.join(parts)
        else:
            self.export_name = self.set_name
    
    def to_dict(self):
        return {
            'set_id': self.set_id,
            'set_name': self.set_name,
            'original_name': self.original_name,
            'spec_name': self.spec_name,
            'diameter_mm': self.diameter_mm,
            'element_ids': self.element_ids,
            'manual': self.manual,
            'problematic': self.problematic,
            'diameter_missing': self.diameter_missing,
            'export_name': self.export_name
        }
   
    @classmethod
    def from_dict(cls, data):
        group = cls(
            data.get('set_id', 0),
            data.get('set_name', ''),
            data.get('spec_name', ''),
            data.get('diameter_mm', 6.35)
        )
        group.original_name = data.get('original_name', data.get('set_name', ''))
        group.element_ids = data.get('element_ids', [])
        group.manual = data.get('manual', False)
        group.problematic = data.get('problematic', False)
        group.diameter_missing = data.get('diameter_missing', False)
        group.export_name = data.get('export_name', data.get('set_name', ''))
        return group

class FastenerGroupDialog(QDialog):
    def __init__(self, group=None, specs=None, parent=None):
        super().__init__(parent)
        self.group = group or FastenerGroup(set_name="Manual Group")
        self.group.manual = True
        self.specs = specs or {}
        
        self.setWindowTitle("Edit Fastener Group")
        self.setMinimumWidth(500)
        
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        name_layout = QHBoxLayout()
        name_label = QLabel("Group Name:")
        self.name_edit = QLineEdit(self.group.set_name)
        name_layout.addWidget(name_label)
        name_layout.addWidget(self.name_edit)
        layout.addLayout(name_layout)
        
        param_group = QGroupBox("Fastener Parameters")
        param_layout = QFormLayout()
        
        self.diameter_spin = QDoubleSpinBox()
        self.diameter_spin.setRange(1.0, 50.0)
        self.diameter_spin.setDecimals(2)
        self.diameter_spin.setValue(self.group.diameter_mm)
        self.diameter_spin.setSuffix(" mm")
        self.diameter_spin.valueChanged.connect(self._on_diameter_changed)
        param_layout.addRow("Fastener Diameter:", self.diameter_spin)
        
        self.spec_combo = QComboBox()
        self.spec_combo.addItem("")
        if self.specs:
            self.spec_combo.addItems(list(self.specs.keys()))
            if self.group.spec_name in self.specs:
                self.spec_combo.setCurrentText(self.group.spec_name)
        self.spec_combo.currentTextChanged.connect(self._on_spec_changed)
        param_layout.addRow("Fastener Specification:", self.spec_combo)
        
        param_group.setLayout(param_layout)
        layout.addWidget(param_group)
        
        elements_group = QGroupBox("Element Assignment")
        elements_layout = QVBoxLayout()
        
        self.elements_edit = QTextEdit()
        self.elements_edit.setMaximumHeight(100)
        self.elements_edit.setText(", ".join(map(str, self.group.element_ids)))
        self.elements_edit.setPlaceholderText("Enter comma-separated or space-separated element IDs")
        elements_layout.addWidget(QLabel("Element IDs:"))
        elements_layout.addWidget(self.elements_edit)
        
        elements_group.setLayout(elements_layout)
        layout.addWidget(elements_group)
        
        button_layout = QHBoxLayout()
        save_button = QPushButton("Save")
        save_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        button_layout.addWidget(save_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
    
    def _on_diameter_changed(self):
        if not self.group.manual:
            self.group.diameter_mm = self.diameter_spin.value()
            self.group.update_set_name()
            self.name_edit.setText(self.group.set_name)
    
    def _on_spec_changed(self):
        if not self.group.manual:
            self.group.spec_name = self.spec_combo.currentText()
            self.group.update_set_name()
            self.name_edit.setText(self.group.set_name)
    
    def get_group(self):
        self.group.set_name = self.name_edit.text()
        self.group.diameter_mm = self.diameter_spin.value()
        self.group.spec_name = self.spec_combo.currentText()
        
        self.group.update_export_name()
        
        element_text = self.elements_edit.toPlainText().strip()
        if element_text:
            elements = []
            for item in element_text.replace(',', ' ').split():
                try:
                    elements.append(int(item))
                except ValueError:
                    pass
            self.group.element_ids = elements
        else:
            self.group.element_ids = []
        
        return self.group

class FastenerSpecDialog(QDialog):
    def __init__(self, spec=None, parent=None):
        super().__init__(parent)
        self.spec = spec or FastenerSpec()
        self.setWindowTitle(f"Edit Fastener Specification: {self.spec.spec_name}")
        self.setMinimumWidth(400)
        
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        name_layout = QHBoxLayout()
        name_label = QLabel("Specification Name:")
        self.name_edit = QLineEdit(self.spec.spec_name)
        self.name_edit.setReadOnly(bool(self.spec.spec_name))
        name_layout.addWidget(name_label)
        name_layout.addWidget(self.name_edit)
        layout.addLayout(name_layout)
        
        param_group = QGroupBox("Fastener Parameters")
        param_layout = QFormLayout()
        
        self.type_combo = QComboBox()
        self.type_combo.addItems(["Bolt", "Rivet", "Hi-Lok", "Blind Rivet"])
        if self.spec.fastener_type in ["Bolt", "Rivet", "Hi-Lok", "Blind Rivet"]:
            self.type_combo.setCurrentText(self.spec.fastener_type)
        param_layout.addRow("Fastener Type:", self.type_combo)
        
        self.material_combo = QComboBox()
        self.material_combo.addItems(["Titanium","Steel", "Aluminum", "Custom"])
        if self.spec.material in ["Titanium", "Steel", "Aluminum", "Custom"]:
            self.material_combo.setCurrentText(self.spec.material)
        self.material_combo.currentTextChanged.connect(self._update_material)
        param_layout.addRow("Fastener Material:", self.material_combo)
        
        self.elasticity_spin = QDoubleSpinBox()
        self.elasticity_spin.setRange(10.0, 500.0)
        self.elasticity_spin.setSingleStep(5.0)
        self.elasticity_spin.setDecimals(1)
        self.elasticity_spin.setValue(self.spec.elastic_modulus)
        self.elasticity_spin.setSuffix(" GPa")
        self.elasticity_spin.setEnabled(self.spec.material == "Custom")
        param_layout.addRow("Elastic Modulus:", self.elasticity_spin)
        
        param_group.setLayout(param_layout)
        layout.addWidget(param_group)
        
        button_layout = QHBoxLayout()
        save_button = QPushButton("Save")
        save_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        button_layout.addWidget(save_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
        
        self._update_material(self.material_combo.currentText())
    
    def _update_material(self, material):
        if material == "Steel":
            self.elasticity_spin.setValue(210.0)
            self.elasticity_spin.setEnabled(False)
        elif material == "Titanium":
            self.elasticity_spin.setValue(110.0)
            self.elasticity_spin.setEnabled(False)
        elif material == "Aluminum":
            self.elasticity_spin.setValue(70.0)
            self.elasticity_spin.setEnabled(False)
        else:
            self.elasticity_spin.setEnabled(True)
    
    def get_spec(self):
        self.spec.spec_name = self.name_edit.text()
        self.spec.fastener_type = self.type_combo.currentText()
        self.spec.material = self.material_combo.currentText()
        self.spec.elastic_modulus = self.elasticity_spin.value()
        return self.spec

class BdfParser:
    def __init__(self, filename):
        self.filename = filename
        self.sets = {}
        self.set_names = {}
        self.fastener_groups = []
    
    def parse(self):
        try:
            with open(self.filename, 'r') as file:
                content = file.read()
            
            clean_content = self._preprocess_content(content)
            lines = clean_content.splitlines()
            
            current_set_id = None
            current_elements = []
            
            for line in lines:
                line = line.strip()
                
                if line.startswith('SET '):
                    if current_set_id is not None:
                        self.sets[current_set_id] = current_elements
                    
                    parts = line.split('=')
                    if len(parts) == 2:
                        set_id_part = parts[0].strip()
                        set_id_match = re.search(r'SET\s+(\d+)', set_id_part)
                        if set_id_match:
                            current_set_id = int(set_id_match.group(1))
                            
                            element_part = parts[1].strip()
                            elements = self._parse_element_list(element_part)
                            current_elements = elements
                
                elif line.startswith('$HMSET '):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            set_id = int(parts[1])
                            set_name_match = re.search(r'"([^"]*)"', line)
                            if set_name_match:
                                set_name = set_name_match.group(1)
                                self.set_names[set_id] = set_name
                        except (ValueError, IndexError):
                            pass
            
            if current_set_id is not None:
                self.sets[current_set_id] = current_elements
            
            self._create_fastener_groups()
            
            return True
            
        except Exception as e:
            print(f"Error parsing BDF file: {str(e)}")
            return False
    
    def _preprocess_content(self, content):
        content = re.sub(r',\s*\n\s+', ', ', content)
        return content
    
    def _parse_element_list(self, element_text):
        elements = []
        for item in element_text.split(','):
            item = item.strip()
            if item:
                try:
                    elements.append(int(item))
                except ValueError:
                    if "THRU" in item:
                        for range_item in range(int(item.split("THRU")[0].strip()),int(item.split("THRU")[1].strip())+1,1):
                            elements.append(int(range_item))
                    else:
                        pass
        return elements
    
    def _create_fastener_groups(self):
        for set_id, elements in self.sets.items():
            if set_id in self.set_names:
                set_name = self.set_names[set_id]
                
                group = FastenerGroup(set_id, set_name)
                group.original_name = set_name
                group.element_ids = elements
                
                self._parse_fastener_info(group, set_name)
                group.update_export_name()
                
                self.fastener_groups.append(group)
        
        return self.fastener_groups
    
    def _parse_fastener_info(self, group, set_name):
        try:
            parts = set_name.split('__')
            
            if len(parts) >= 3:
                group.diameter_mm = float(parts[2])
                
                if len(parts) >= 4:
                    try:
                        group.spec_name = parts[3]
                        
                    except ValueError:
                        group.problematic = True
                        group.spec_name = "Default"
                else:
                    group.diameter_missing = True
                
                group.update_set_name()
                        
            else:
                group.problematic = True
                group.diameter_missing = True
                
        except Exception:
            group.problematic = True
            group.diameter_missing = True
    
    def export_sets_bdf(self, filename, fastener_groups):
        try:
            with open(filename, 'w') as f:
                f.write("$ Exported Fastener Sets\n")
                f.write("$ Generated by CBUSH Huth Stiffness Updater\n")
                f.write("$\n")
                
                for group in fastener_groups:
                    if group.element_ids:
                        set_id_str = f"{group.set_id}"
                        
                        f.write(f"SET {set_id_str:<8} = ")
                        
                        elements_per_line = 6
                        for i, elem_id in enumerate(group.element_ids):
                            if i > 0 and i % elements_per_line == 0:
                                f.write(",\n        ")
                            elif i > 0:
                                f.write(",")
                            f.write(str(elem_id))
                        f.write("\n")
                        
                        f.write(f"$HMSET {set_id_str:>8} {2:>8} \"{group.export_name}\" 18\n")
                        f.write(f"$HMSETTYPE{set_id_str:>8} \"regular\" 18\n")
                        f.write("$\n")
                
            return True
        except Exception as e:
            print(f"Error exporting sets BDF: {str(e)}")
            return False

class ParallelCalculationWorker(QThread):
    progress_updated = pyqtSignal(int)
    calculation_finished = pyqtSignal(bool, str)
    log_message = pyqtSignal(str)
    groups_updated = pyqtSignal()  # New signal to update GUI
    
    def __init__(self, input_file, output_file, default_parameters, fastener_groups, 
                 specs_dict, other_elements, num_workers, connection_mode):
        super().__init__()
        self.input_file = input_file
        self.output_file = output_file
        self.default_parameters = default_parameters
        self.fastener_groups = fastener_groups
        self.specs_dict = specs_dict
        self.other_elements = other_elements
        self.num_workers = num_workers
        self.connection_mode = connection_mode  # 'convert_to_airbus', 'convert_to_normal', 'mix_mode'
    
    def _validate_element_assignments(self, elements_to_update, connections):
        """Validate that all elements are properly assigned to connections"""
        element_to_connection = {}
        for conn_id, connection in connections.items():
            for elem_id in connection.cbush_elements:
                element_to_connection[elem_id] = conn_id
        
        unassigned = []
        for elem_id in elements_to_update:
            if elem_id not in element_to_connection:
                unassigned.append(elem_id)
        
        if unassigned:
            self.log_message.emit(f"ERROR: {len(unassigned)} elements not assigned to connections: {unassigned}")
            return False
        
        self.log_message.emit(f"Validation passed: All {len(elements_to_update)} elements assigned to connections")
        return True
    
    def _validate_group_assignments(self, connections):
        """Validate and fix group assignments to ensure complete connections"""
        self.log_message.emit("Validating and updating fastener group assignments...")
        
        # Create element to connection mapping
        element_to_connection = {}
        for conn_id, connection in connections.items():
            for elem_id in connection.cbush_elements:
                element_to_connection[elem_id] = conn_id
        
        # Create element to group mapping
        element_to_group = {}
        for group in self.fastener_groups:
            for elem_id in group.element_ids:
                if elem_id in element_to_group:
                    self.log_message.emit(f"ERROR: Element {elem_id} is in multiple groups!")
                    return False
                element_to_group[elem_id] = group
        
        # Find connections that are partially in groups
        connection_groups = {}  # conn_id -> set of groups
        for elem_id, conn_id in element_to_connection.items():
            if elem_id in element_to_group:
                group = element_to_group[elem_id]
                if conn_id not in connection_groups:
                    connection_groups[conn_id] = set()
                connection_groups[conn_id].add(group)
        
        # Check for connections in multiple groups
        groups_updated = False
        for conn_id, groups in connection_groups.items():
            if len(groups) > 1:
                self.log_message.emit(f"ERROR: Connection {conn_id} spans multiple groups: {[g.set_name for g in groups]}")
                return False
            
            # Check if connection is complete in the group
            connection = connections[conn_id]
            group = list(groups)[0]
            missing_elements = []
            for elem_id in connection.cbush_elements:
                if elem_id not in group.element_ids:
                    missing_elements.append(elem_id)
            
            if missing_elements:
                self.log_message.emit(f"Adding {len(missing_elements)} missing elements to group {group.set_name} for complete connection {conn_id}")
                group.element_ids.extend(missing_elements)
                groups_updated = True
        
        if groups_updated:
            self.groups_updated.emit()
            self.log_message.emit("Updated fastener groups to include complete connections")
        
        return True
    
    def _convert_airbus_to_normal(self, model, connections, model_data):
        """Convert Airbus connections to normal connections"""
        self.log_message.emit("Converting Airbus connections to normal connections...")
        
        converted_connections = 0
        elements_to_remove = []
        properties_to_remove = []
        
        for connection in connections.values():
            if not connection.is_airbus_method or not connection.airbus_cbushes:
                continue
            
            # Remove Airbus CBUSHes and their properties
            for airbus_cbush_id in connection.airbus_cbushes:
                if airbus_cbush_id in model.elements:
                    airbus_elem = model.elements[airbus_cbush_id]
                    prop_id = airbus_elem.pid
                    
                    # Remove element
                    del model.elements[airbus_cbush_id]
                    elements_to_remove.append(airbus_cbush_id)
                    
                    # Remove property if not used elsewhere
                    prop_used_elsewhere = any(
                        elem.pid == prop_id for elem_id, elem in model.elements.items() 
                        if elem_id != airbus_cbush_id and hasattr(elem, 'pid')
                    )
                    
                    if not prop_used_elsewhere and prop_id in model.properties:
                        del model.properties[prop_id]
                        properties_to_remove.append(prop_id)
            
            # Update connection
            connection.airbus_cbushes = []
            connection.airbus_cbush = None
            connection.is_airbus_method = False
            connection.cbush_elements = [elem_id for elem_id in connection.cbush_elements 
                                        if elem_id not in elements_to_remove]
        
        converted_connections += 1
        
        # Remove elements from fastener groups
        if elements_to_remove:
            for group in self.fastener_groups:
                original_count = len(group.element_ids)
                group.element_ids = [elem_id for elem_id in group.element_ids 
                                    if elem_id not in elements_to_remove]
                removed_count = original_count - len(group.element_ids)
                if removed_count > 0:
                    self.log_message.emit(f"Removed {removed_count} Airbus elements from group {group.set_name}")
        
        self.log_message.emit(f"Converted {converted_connections} Airbus connections to normal")
        self.log_message.emit(f"Removed {len(elements_to_remove)} Airbus CBUSH elements")
        self.log_message.emit(f"Removed {len(properties_to_remove)} unused PBUSH properties")
        
        if elements_to_remove:
            self.groups_updated.emit()
        
        return elements_to_remove
    
    def run(self):
        try:
            from pyNastran.bdf.bdf import BDF
            
            self.log_message.emit("Loading BDF file...")
            self.progress_updated.emit(5)
            
            model = try_bdf_file(self.input_file)
            
            self.log_message.emit("Identifying CBUSH elements...")
            self.progress_updated.emit(10)
            
            cbush_elements = {}
            for elem_id, elem in model.elements.items():
                if elem.type == 'CBUSH':
                    cbush_elements[elem_id] = elem
            
            self.log_message.emit(f"Found {len(cbush_elements)} CBUSH elements.")
            
            model_data = self._prepare_model_data(model)
            
            # Determine which elements to update
            elements_to_update = set()
            
            for group in self.fastener_groups:
                elements_to_update.update(group.element_ids)
            
            if self.other_elements:
                elements_to_update.update(self.other_elements)
            elif not self.fastener_groups:
                elements_to_update = set(cbush_elements.keys())
            
            elements_to_update = elements_to_update.intersection(set(cbush_elements.keys()))
            
            self.log_message.emit(f"Will update {len(elements_to_update)} CBUSH elements.")
            
            if not elements_to_update:
                self.log_message.emit("No CBUSH elements to update.")
                self.calculation_finished.emit(False, "No elements to update.")
                return
            
            self.log_message.emit("Identifying multi-plate connections in parallel...")
            self.progress_updated.emit(15)
            
            connections = self._parallel_identify_connections(list(elements_to_update), model_data)
            
            # Store model_data reference in connections for material orientation calculations
            for connection in connections.values():
                connection._model_data = model_data
            
            # Validate and update group assignments
            if not self._validate_group_assignments(connections):
                self.calculation_finished.emit(False, "Group assignment validation failed")
                return
            
            if not self._validate_element_assignments(elements_to_update, connections):
                self.calculation_finished.emit(False, "Element assignment validation failed")
                return
            
            self.log_message.emit(f"Found {len(connections)} multi-plate connections.")
            
            self.log_message.emit("Extracting plate properties...")
            self.progress_updated.emit(25)
            
            analyzer = CBUSHAnalyzer(model)
            for connection in connections.values():
                for plate_id in connection.plates:
                    if plate_id not in connection.plate_properties:
                        # Use simple plate properties extraction WITHOUT orientation transformation
                        connection.plate_properties[plate_id] = analyzer._get_simple_plate_properties(plate_id)
                connection.finalize_with_airbus_detection(model_data)
            
            # Handle conversion modes
            new_cbush_elements = {}
            
            if self.connection_mode == 'convert_to_normal':
                self.log_message.emit("Converting Airbus connections to normal...")
                self.progress_updated.emit(30)
                removed_elements = self._convert_airbus_to_normal(model, connections, model_data)
                # Update elements_to_update to exclude removed elements
                elements_to_update = elements_to_update - set(removed_elements)
                
            elif self.connection_mode == 'convert_to_airbus':
                self.log_message.emit("Processing Airbus method connections...")
                self.progress_updated.emit(30)
                new_cbush_elements = self._process_airbus_method(model, connections, model_data)
                self._update_fastener_groups_with_new_elements(new_cbush_elements, connections)
                
            else:  # mix_mode
                self.log_message.emit("Using mix mode - preserving existing connection types...")
                self.progress_updated.emit(30)
            
            element_to_group = {}
            for group in self.fastener_groups:
                for elem_id in group.element_ids:
                    element_to_group[elem_id] = group
            
            processed_pids = set()
            next_pid = max(model.properties.keys()) + 1 if model.properties else 1
            
            self.log_message.emit("Preparing element parameters...")
            self.progress_updated.emit(35)
            
            element_params = []
            config_usage = {group.set_name: 0 for group in self.fastener_groups}
            config_usage["Default"] = 0
            
            element_to_connection = {}
            for conn_id, connection in connections.items():
                for elem_id in connection.cbush_elements:
                    element_to_connection[elem_id] = conn_id
            
            all_elements_to_update = elements_to_update.copy()
            for new_elem_id in new_cbush_elements.keys():
                all_elements_to_update.add(new_elem_id)
            
            for elem_id in all_elements_to_update:
                connection_id = element_to_connection.get(elem_id)
                
                if elem_id in element_to_group:
                    group = element_to_group[elem_id]
                    
                    if group.spec_name in self.specs_dict:
                        spec = self.specs_dict[group.spec_name]
                        fastener_type = spec.fastener_type.lower()
                        fastener_e = spec.elastic_modulus
                    else:
                        fastener_type = self.default_parameters['fastener_type'].lower()
                        fastener_e = self.default_parameters['fastener_e']
                    
                    # Determine use_airbus_method based on connection mode
                    if self.connection_mode == 'convert_to_airbus':
                        use_airbus_method = True
                    elif self.connection_mode == 'convert_to_normal':
                        use_airbus_method = False
                    else:  # mix_mode
                        # Use existing connection type
                        connection = connections.get(connection_id)
                        use_airbus_method = connection.is_airbus_method if connection else False
                    
                    parameters = {
                        'fastener_diameter_mm': group.diameter_mm,
                        'fastener_type': fastener_type,
                        'fastener_e': fastener_e,
                        'connection_method': self.default_parameters['connection_method'],
                        'k4': self.default_parameters['k4'],
                        'k5': self.default_parameters['k5'],
                        'k6': self.default_parameters['k6'],
                        'use_airbus_method': use_airbus_method,
                        'use_material_orientation': self.default_parameters.get('use_material_orientation', False)
                    }
                    config_usage[group.set_name] += 1
                else:
                    # Determine use_airbus_method for default parameters
                    if self.connection_mode == 'convert_to_airbus':
                        use_airbus_method = True
                    elif self.connection_mode == 'convert_to_normal':
                        use_airbus_method = False
                    else:  # mix_mode
                        connection = connections.get(connection_id)
                        use_airbus_method = connection.is_airbus_method if connection else False
                    
                    parameters = self.default_parameters.copy()
                    parameters['use_airbus_method'] = use_airbus_method
                    config_usage["Default"] += 1
                
                element_params.append((elem_id, connection_id, parameters))
            
            self.log_message.emit("Calculating stiffness values in parallel...")
            self.progress_updated.emit(50)
            
            stiffness_results = self._parallel_calculate_stiffness(
                element_params, model_data, connections, self.specs_dict)
            
            self.log_message.emit("Updating model with calculated stiffness...")
            self.progress_updated.emit(70)
            
            for elem_id in all_elements_to_update:
                if elem_id not in model.elements:
                    continue
                    
                elem = model.elements[elem_id]
                original_pid = elem.pid
                
                if original_pid in processed_pids:
                    new_pid = next_pid
                    next_pid += 1
                    original_prop = model.properties[original_pid]
                    from pyNastran.bdf.cards.properties.bush import PBUSH
                    new_prop = PBUSH(new_pid, original_prop.Ki.copy() if hasattr(original_prop, 'Ki') else [0]*6, original_prop.Bi.copy() if hasattr(original_prop, 'Bi') else [0]*6, original_prop.GEi.copy() if hasattr(original_prop, 'GEi') else [0]*6, original_prop.sa if hasattr(original_prop, 'sa') else None, original_prop.st if hasattr(original_prop, 'st') else None)
                    model.properties[new_pid] = new_prop
                    elem.pid = new_pid
                    current_pid = new_pid
                else:
                    current_pid = original_pid
                
                processed_pids.add(current_pid)
                
                if elem_id in stiffness_results and current_pid in model.properties:
                    pbush = model.properties[current_pid]
                    stiffness = stiffness_results[elem_id]
                    
                    for k, v in stiffness.items():
                        if k == 'K1':
                            pbush.Ki[0] = v
                        elif k == 'K2':
                            pbush.Ki[1] = v
                        elif k == 'K3':
                            pbush.Ki[2] = v
                        elif k == 'K4':
                            pbush.Ki[3] = v
                        elif k == 'K5':
                            pbush.Ki[4] = v
                        elif k == 'K6':
                            pbush.Ki[5] = v
            
            self.log_message.emit("\nElement configuration summary:")
            for config_name, count in config_usage.items():
                if count > 0:
                    self.log_message.emit(f"  - {config_name}: {count} elements")
            
            if self.connection_mode == 'convert_to_airbus' and new_cbush_elements:
                self.log_message.emit(f"\nAirbus method: Created {len(new_cbush_elements)} new CBUSH elements")
                airbus_detected = sum(1 for conn in connections.values() if conn.is_airbus_method)
                self.log_message.emit(f"Airbus method: Detected {airbus_detected} existing Airbus connections")
            elif self.connection_mode == 'mix_mode':
                airbus_connections = sum(1 for conn in connections.values() if conn.is_airbus_method)
                normal_connections = len(connections) - airbus_connections
                self.log_message.emit(f"\nMix mode: {airbus_connections} Airbus connections, {normal_connections} normal connections")
            
            self.log_message.emit(f"\nSaving modified BDF to {self.output_file}...")
            self.progress_updated.emit(85)
            model.write_bdf(self.output_file)
            
            self.progress_updated.emit(100)
            self.log_message.emit("Parallel calculation completed successfully!")
            
            # Update fastener groups table after completion if new elements were created
            if self.connection_mode in ['convert_to_airbus', 'convert_to_normal']:
                self.groups_updated.emit()
            
            if self.connection_mode == 'convert_to_airbus' and new_cbush_elements:
                self.log_message.emit("Note: Fastener groups have been updated with new Airbus CBUSH elements. Export updated sets to include them.")
            elif self.connection_mode == 'convert_to_normal':
                self.log_message.emit("Note: Airbus CBUSH elements have been removed from fastener groups. Export updated sets to reflect changes.")
            
            self.calculation_finished.emit(True, "CBUSH stiffness values updated successfully!")
            
        except Exception as e:
            self.log_message.emit(f"Error: {str(e)}")
            self.calculation_finished.emit(False, f"Error: {str(e)}")
    
    def _prepare_model_data(self, model):
        model_data = {
            'elements': dict(model.elements),
            'rigid_elements': dict(model.rigid_elements),
            'properties': dict(model.properties),
            'materials': dict(model.materials),
            'nodes': dict(model.nodes),  # Added nodes for coordinate calculations
            'coords': dict(model.coords) if hasattr(model, 'coords') else {}  # Added coordinate systems
        }
        return model_data
    
    def _parallel_identify_connections(self, cbush_element_ids, model_data):
        chunk_size = max(1, len(cbush_element_ids) // self.num_workers)
        chunks = [cbush_element_ids[i:i + chunk_size] 
                  for i in range(0, len(cbush_element_ids), chunk_size)]
        
        args_list = []
        for i, chunk in enumerate(chunks):
            args_list.append((chunk, model_data))
        
        all_connections = {}
        all_processed = set()
        
        with ProcessPoolExecutor(max_workers=self.num_workers) as executor:
            futures = [executor.submit(process_cbush_connection_chunk, args) for args in args_list]
            
            for future in as_completed(futures):
                try:
                    connections, processed = future.result()
                    for conn_id, connection in connections.items():
                        if not any(elem_id in all_processed for elem_id in connection.cbush_elements):
                            new_conn_id = len(all_connections) + 1
                            connection.connection_id = new_conn_id
                            all_connections[new_conn_id] = connection
                            all_processed.update(connection.cbush_elements)
                except Exception as e:
                    self.log_message.emit(f"Error in connection identification: {str(e)}")
        
        return all_connections
    
    def _parallel_calculate_stiffness(self, element_params, model_data, connections, specs_dict):
        chunk_size = max(1, len(element_params) // self.num_workers)
        chunks = [element_params[i:i + chunk_size] 
                  for i in range(0, len(element_params), chunk_size)]
        
        args_list = []
        for chunk in chunks:
            args_list.append((chunk, model_data, connections, self.default_parameters, specs_dict))
        
        all_stiffness_results = {}
        
        with ProcessPoolExecutor(max_workers=self.num_workers) as executor:
            futures = [executor.submit(calculate_stiffness_chunk, args) for args in args_list]
            
            completed_chunks = 0
            for future in as_completed(futures):
                try:
                    stiffness_results = future.result()
                    all_stiffness_results.update(stiffness_results)
                    completed_chunks += 1
                    progress = 50 + (15 * completed_chunks // len(chunks))
                    self.progress_updated.emit(progress)
                except Exception as e:
                    self.log_message.emit(f"Error in stiffness calculation: {str(e)}")
        
        return all_stiffness_results
    
    def _process_airbus_method(self, model, connections, model_data):
        from pyNastran.bdf.cards.elements.bush import CBUSH
        from pyNastran.bdf.cards.properties.bush import PBUSH
        
        new_cbush_elements = {}
        next_elem_id = max(model.elements.keys()) + 1 if model.elements else 1
        next_prop_id = max(model.properties.keys()) + 1 if model.properties else 1
        
        converted_connections = 0
        
        for connection in connections.values():
            if len(connection.plates) <= 2:
                continue
            
            if connection.is_airbus_method:
                self.log_message.emit(f"Connection {connection.connection_id} already uses Airbus method - skipping")
                continue
            
            outer_nodes = connection.get_outer_nodes_for_new_airbus(model_data)
            if len(outer_nodes) != 2:
                self.log_message.emit(f"Could not determine outer nodes for connection {connection.connection_id}")
                continue
            
            new_elem_id = next_elem_id
            next_elem_id += 1
            
            new_prop_id = next_prop_id
            next_prop_id += 1
            
            new_pbush = PBUSH(new_prop_id, [1e6, 100, 100, 100, 1e8, 1e8], None, None)
            model.properties[new_prop_id] = new_pbush
            
            new_cbush = CBUSH(new_elem_id, new_prop_id, outer_nodes, x=[1.0,0.0,0.0], g0=None)
            model.elements[new_elem_id] = new_cbush
            
            # Update connection to include the new Airbus CBUSH
            connection.cbush_elements.append(new_elem_id)
            connection.airbus_cbushes.append(new_elem_id)
            connection.airbus_cbush = new_elem_id
            connection.is_airbus_method = True
            connection.end_nodes = outer_nodes
            
            new_cbush_elements[new_elem_id] = connection.connection_id
            converted_connections += 1
            
            self.log_message.emit(f"Created Airbus CBUSH {new_elem_id} for connection {connection.connection_id}")
        
        self.log_message.emit(f"Converted {converted_connections} connections to Airbus method")
        return new_cbush_elements
    
    def _update_fastener_groups_with_new_elements(self, new_cbush_elements, connections):
        connection_to_groups = {}
        
        for group in self.fastener_groups:
            for elem_id in group.element_ids:
                for conn_id, connection in connections.items():
                    if elem_id in connection.regular_cbushes:
                        if conn_id not in connection_to_groups:
                            connection_to_groups[conn_id] = []
                        if group not in connection_to_groups[conn_id]:
                            connection_to_groups[conn_id].append(group)
                        break
        
        for new_elem_id, conn_id in new_cbush_elements.items():
            if conn_id in connection_to_groups:
                for group in connection_to_groups[conn_id]:
                    if new_elem_id not in group.element_ids:
                        group.element_ids.append(new_elem_id)
                        self.log_message.emit(f"Added Airbus CBUSH {new_elem_id} to group {group.set_name}")

class CBUSHAnalyzer:
    def __init__(self, model):
        self.model = model
        self.node_to_cbush = defaultdict(list)
        self.cbush_to_plates = {}
        self.plate_properties = {}
        
        for elem_id, elem in model.elements.items():
            if elem.type == 'CBUSH':
                for node_id in elem.nodes:
                    self.node_to_cbush[node_id].append(elem_id)
    
    def _get_cbush_plates(self, cbush_id):
        cbush_elem = self.model.elements[cbush_id]
        node1, node2 = cbush_elem.nodes
        
        plates_node1 = self._find_connected_plates(node1)
        plates_node2 = self._find_connected_plates(node2)
        
        plates = []
        if plates_node1:
            plates.append(plates_node1[0])
        if plates_node2:
            plates.append(plates_node2[0])
        
        return plates
    
    def _find_connected_plates(self, node_id):
        connected_plates = []
        
        rbe3_independent_nodes = []
        for elem_id, elem in self.model.rigid_elements.items():
            if hasattr(elem, 'refgrid') and elem.refgrid == node_id:
                if hasattr(elem, 'Gijs') and len(elem.Gijs) > 0:
                    if len(elem.Gijs[0]) > 0:
                        rbe3_independent_nodes.append(elem.Gijs[0][0])
                break
        
        if rbe3_independent_nodes:
            plates = self._find_direct_connected_plates(rbe3_independent_nodes[0])
            connected_plates.extend(plates)
        else:
            plates = self._find_direct_connected_plates(node_id)
            connected_plates.extend(plates)
        
        return connected_plates
    
    def _find_direct_connected_plates(self, node_id):
        connected_plates = []
        
        for elem_id, elem in self.model.elements.items():
            if elem.type in ['CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6']:
                if node_id in elem.node_ids:
                    connected_plates.append(elem_id)
                    break
        
        return connected_plates
    
    def _get_simple_plate_properties(self, elem_id):
        """Enhanced plate properties with full composite support but NO orientation transformation"""
        plate_props = {
            'type': 'metal',
            'thickness': 1.0,
            'elastic_modulus': 70000.0,
            'elastic_modulus_x': 70000.0,
            'elastic_modulus_y': 70000.0,
            'shear_modulus_xy': 26900.0,
            'use_directional': False
        }
        
        if elem_id in self.model.elements:
            elem = self.model.elements[elem_id]
            
            if elem.type in ['CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6']:
                prop_id = elem.pid
                
                if prop_id in self.model.properties:
                    prop = self.model.properties[prop_id]
                    
                    if prop.type == 'PSHELL':
                        return self._handle_pshell_properties(prop, elem)
                    
                    elif prop.type in ['PCOMP', 'PCOMPG', 'PCOMPP']:
                        return self._handle_composite_properties(prop, elem)
        
        return plate_props
    
    def _handle_pshell_properties(self, prop, elem):
        """Handle PSHELL properties without orientation transformation"""
        thickness = prop.t
        mid1 = prop.mid1
        
        # Get basic material properties without transformation
        material_props = self._get_basic_material_properties(mid1)
        
        return {
            'type': material_props.get('type', 'metal'),
            'id': elem.eid,
            'property_id': prop.pid,
            'thickness': thickness,
            'material_id': mid1,
            'elastic_modulus': material_props.get('E_effective', 70000.0),
            'elastic_modulus_x': material_props.get('Ex', material_props.get('E_effective', 70000.0)),
            'elastic_modulus_y': material_props.get('Ey', material_props.get('E_effective', 70000.0)),
            'shear_modulus_xy': material_props.get('Gxy', 26900.0),
            'use_directional': False,
            'E1': material_props.get('E1', material_props.get('E_effective', 70000.0)),
            'E2': material_props.get('E2', material_props.get('E_effective', 70000.0)),
            'G12': material_props.get('G12', 26900.0),
            'NU12': material_props.get('NU12', 0.3)
        }
    
    def _handle_composite_properties(self, prop, elem):
        """Handle composite properties (PCOMP, PCOMPG, PCOMPP) with full CLT"""
        
        # Extract plies based on property type
        if prop.type == 'PCOMP':
            plies = self._extract_pcomp_plies(prop)
        elif prop.type == 'PCOMPG':
            plies = self._extract_pcompg_plies(prop, elem)
        elif prop.type == 'PCOMPP':
            plies = self._extract_pcompp_plies(prop, elem)
        else:
            plies = []
        
        if not plies:
            return self._get_default_plate_props()
        
        # Always use CLT for composites (single or multi-ply)
        effective_props = self._calculate_clt_effective_properties(plies)
        
        return {
            'type': 'composite',
            'id': elem.eid,
            'property_id': prop.pid,
            'thickness': effective_props['total_thickness'],
            'plies': plies,
            'elastic_modulus': effective_props['E_effective'],
            'elastic_modulus_x': effective_props['Ex'],
            'elastic_modulus_y': effective_props['Ey'],
            'shear_modulus_xy': effective_props['Gxy'],
            'use_directional': False,  # No transformation here
            'clt_matrices': effective_props.get('clt_matrices', {}),
            'E1': effective_props['Ex'],
            'E2': effective_props['Ey'],
            'G12': effective_props['Gxy'],
            'NU12': effective_props.get('NU12', 0.3)
        }
    
    def _extract_pcomp_plies(self, prop):
        """Extract plies from PCOMP property"""
        plies = []
        
        try:
            n_plies = len(prop.mids) if hasattr(prop, 'mids') else 0
            
            # Handle symmetry
            is_symmetric = hasattr(prop, 'lam') and prop.lam in ['SYM', 'SYMEM', 'SYBEND', 'SYSMEAR']
            
            # Basic plies
            for i in range(n_plies):
                if i < len(prop.mids) and prop.mids[i] is not None:
                    ply = {
                        'material_id': prop.mids[i],
                        'thickness': prop.thicknesses[i] if i < len(prop.thicknesses) else 0.1,
                        'angle': prop.thetas[i] if i < len(prop.thetas) else 0.0,
                        'ply_id': i + 1
                    }
                    plies.append(ply)
            
            # Add symmetric plies if needed
            if is_symmetric:
                sym_plies = []
                for ply in reversed(plies):
                    sym_ply = ply.copy()
                    sym_ply['ply_id'] = len(plies) + len(sym_plies) + 1
                    sym_plies.append(sym_ply)
                plies.extend(sym_plies)
                
        except Exception as e:
            print(f"Error extracting PCOMP plies: {e}")
        
        return plies
    
    def _extract_pcompg_plies(self, prop, elem):
        """Extract plies from PCOMPG property (global plies)"""
        plies = []
        
        try:
            # PCOMPG references global plies
            if hasattr(prop, 'global_ply_ids'):
                for ply_id in prop.global_ply_ids:
                    if ply_id in self.model.plies:
                        global_ply = self.model.plies[ply_id]
                        ply = {
                            'material_id': global_ply.mid,
                            'thickness': global_ply.t,
                            'angle': global_ply.theta,
                            'ply_id': ply_id,
                            'global_ply': True
                        }
                        plies.append(ply)
            
            # Alternative: check element's ply list
            elif hasattr(elem, 'plylist') and elem.plylist:
                for ply_id in elem.plylist:
                    if ply_id in self.model.plies:
                        global_ply = self.model.plies[ply_id]
                        ply = {
                            'material_id': global_ply.mid,
                            'thickness': global_ply.t,
                            'angle': global_ply.theta,
                            'ply_id': ply_id,
                            'global_ply': True
                        }
                        plies.append(ply)
                        
        except Exception as e:
            print(f"Error extracting PCOMPG plies: {e}")
        
        return plies
    
    def _extract_pcompp_plies(self, prop, elem):
        """Extract plies from PCOMPP property (parametric)"""
        plies = []
        
        try:
            # PCOMPP is parametric - implementation depends on specific usage
            if hasattr(prop, 'pcomp_id') and prop.pcomp_id in self.model.properties:
                # Reference to base PCOMP
                base_pcomp = self.model.properties[prop.pcomp_id]
                plies = self._extract_pcomp_plies(base_pcomp)
                
                # Apply parametric modifications if available
                if hasattr(prop, 'thickness_scale'):
                    for ply in plies:
                        ply['thickness'] *= prop.thickness_scale
                        
        except Exception as e:
            print(f"Error extracting PCOMPP plies: {e}")
        
        return plies
    
    def _calculate_clt_effective_properties(self, plies):
        """Calculate effective properties using Classical Lamination Theory"""
        if not plies:
            return {
                'E_effective': 70000.0,
                'Ex': 70000.0,
                'Ey': 70000.0,
                'Gxy': 26900.0,
                'NU12': 0.3,
                'total_thickness': 1.0
            }
        
        total_thickness = sum(ply['thickness'] for ply in plies)
        
        # Initialize CLT matrices
        A = np.zeros((3, 3))  # Extensional stiffness
        B = np.zeros((3, 3))  # Coupling stiffness
        D = np.zeros((3, 3))  # Bending stiffness
        
        z = -total_thickness / 2.0  # Start from bottom
        
        for ply in plies:
            t_ply = ply['thickness']
            mat_id = ply['material_id']
            ply_angle = ply['angle']
            
            # Get material properties
            if mat_id in self.model.materials:
                material = self.model.materials[mat_id]
                
                if material.type == 'MAT8':
                    E1, E2, G12, nu12 = material.e11, material.e22, material.g12, material.nu12
                elif material.type == 'MAT1':
                    E1 = E2 = material.e
                    G12 = material.g
                    nu12 = material.nu
                else:
                    continue
            else:
                continue
            
            # Calculate Q matrix for this ply (ply angle only, no element orientation)
            Q = self._calculate_q_matrix(E1, E2, G12, nu12, np.radians(ply_angle))
            
            z_bottom = z
            z_top = z + t_ply
            
            # Add to ABD matrices
            A += Q * t_ply
            B += Q * (z_top**2 - z_bottom**2) / 2.0
            D += Q * (z_top**3 - z_bottom**3) / 3.0
            
            z = z_top
        
        # Calculate effective properties from A matrix
        # try:
        # Invert A matrix to get compliance
        a = np.linalg.inv(A)
        
        Ex = 1.0 / (a[0, 0] * total_thickness) if a[0, 0] > 0 else 70000.0
        Ey = 1.0 / (a[1, 1] * total_thickness) if a[1, 1] > 0 else 70000.0
        Gxy = 1.0 / (a[2, 2] * total_thickness) if a[2, 2] > 0 else 26900.0
        
        # Calculate effective Poisson's ratio
        nu_xy = -a[0, 1] / a[0, 0] if a[0, 0] > 0 else 0.3
        
        E_effective = (Ex * Ey)**0.5
            
        # except (np.linalg.LinAlgError, ZeroDivisionError):
        #     # Fallback to weighted average
        #     weighted_ex = weighted_ey = weighted_g = weighted_nu = 0.0
        #     total_weight = 0.0
            
        #     for ply in plies:
        #         weight = ply['thickness']
        #         mat_id = ply['material_id']
                
        #         if mat_id in self.model.materials:
        #             material = self.model.materials[mat_id]
        #             if material.type == 'MAT8':
        #                 # Apply simple Jones transformation for weighted average
        #                 ply_angle_rad = np.radians(ply['angle'])
        #                 Ex_ply, Ey_ply, Gxy_ply = self._jones_transformation(
        #                     material.e11, material.e22, material.g12, material.nu12, ply_angle_rad)
        #                 weighted_ex += Ex_ply * weight
        #                 weighted_ey += Ey_ply * weight
        #                 weighted_g += Gxy_ply * weight
        #                 weighted_nu += material.nu12 * weight
        #             elif material.type == 'MAT1':
        #                 weighted_ex += material.e * weight
        #                 weighted_ey += material.e * weight
        #                 weighted_g += material.g * weight
        #                 weighted_nu += material.nu * weight
        #             total_weight += weight
            
        #     if total_weight > 0:
        #         Ex = weighted_ex / total_weight
        #         Ey = weighted_ey / total_weight
        #         Gxy = weighted_g / total_weight
        #         nu_xy = weighted_nu / total_weight
        #         E_effective = (Ex * Ey)**0.5
        #     else:
        #         Ex = Ey = E_effective = 70000.0
        #         Gxy = 26900.0
        #         nu_xy = 0.3
        
        return {
            'E_effective': E_effective,
            'Ex': Ex,
            'Ey': Ey,
            'Gxy': Gxy,
            'NU12': nu_xy,
            'total_thickness': total_thickness,
            'clt_matrices': {'A': A, 'B': B, 'D': D}
        }
    
    def _calculate_q_matrix(self, E1, E2, G12, nu12, theta_rad):
        """Calculate transformed stiffness matrix Q for a ply"""
        # try:
        # Calculate nu21
        nu21 = nu12 * E2 / E1
        
        # Reduced stiffness matrix components in material coordinates
        Q11 = E1 / (1 - nu12 * nu21)
        Q22 = E2 / (1 - nu12 * nu21)
        Q12 = nu12 * E2 / (1 - nu12 * nu21)
        Q66 = G12
        Q16 = Q26 = 0  # No coupling for orthotropic materials
        
        # Transformation parameters
        c = np.cos(theta_rad)
        s = np.sin(theta_rad)
        c2 = c * c
        s2 = s * s
        c4 = c2 * c2
        s4 = s2 * s2
        
        # Transformed stiffness matrix components
        Q = np.zeros((3, 3))
        
        Q[0, 0] = Q11 * c4 + 2 * (Q12 + 2 * Q66) * c2 * s2 + Q22 * s4
        Q[0, 1] = (Q11 + Q22 - 4 * Q66) * c2 * s2 + Q12 * (c4 + s4)
        Q[0, 2] = (Q11 - Q12 - 2 * Q66) * c**3 * s + (Q12 - Q22 + 2 * Q66) * c * s**3
        Q[1, 1] = Q11 * s4 + 2 * (Q12 + 2 * Q66) * c2 * s2 + Q22 * c4
        Q[1, 2] = (Q11 - Q12 - 2 * Q66) * c * s**3 + (Q12 - Q22 + 2 * Q66) * c**3 * s
        Q[2, 2] = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * c2 * s2 + Q66 * (c4 + s4)
        
        # Symmetric components
        Q[1, 0] = Q[0, 1]
        Q[2, 0] = Q[0, 2]
        Q[2, 1] = Q[1, 2]
        
        return Q
    
    def _get_basic_material_properties(self, material_id):
        """Get basic material properties without transformation"""
        material_props = {
            'type': 'default',
            'E_effective': 70000.0,
            'Ex': 70000.0,
            'Ey': 70000.0,
            'use_directional': False
        }
        
        if material_id not in self.model.materials:
            return material_props
        
        material = self.model.materials[material_id]
        
        if material.type == 'MAT1':
            # Isotropic material
            material_props = {
                'type': 'isotropic',
                'E': material.e,
                'G': material.g,
                'nu': material.nu,
                'E_effective': material.e,
                'Ex': material.e,
                'Ey': material.e,
                'use_directional': False,
                'E1': material.e,
                'E2': material.e,
                'G12': material.g,
                'NU12': material.nu
            }
            
        elif material.type == 'MAT8':
            # Orthotropic - no transformation here
            material_props = {
                'type': 'orthotropic',
                'E1': material.e11,
                'E2': material.e22,
                'G12': material.g12,
                'NU12': material.nu12,
                'E_effective': (material.e11 * material.e22)**0.5,
                'Ex': material.e11,
                'Ey': material.e22,
                'use_directional': False
            }
        
        return material_props
    
    def _get_default_plate_props(self):
        """Get default plate properties"""
        return {
            'type': 'metal',
            'thickness': 1.0,
            'elastic_modulus': 70000.0,
            'elastic_modulus_x': 70000.0,
            'elastic_modulus_y': 70000.0,
            'use_directional': False
        }

class HuthCalculatorGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("CBUSH Huth Stiffness Updater")
        self.setMinimumSize(1100, 950)
        
        try:
            if getattr(sys, 'frozen', False):
                app_dir = sys._MEIPASS
            else:
                app_dir = os.path.dirname(os.path.abspath(__file__))
            
            icon_path = os.path.join(app_dir, "huth.ico")
            
            if os.path.exists(icon_path):
                self.setWindowIcon(QIcon(icon_path))
        except Exception as e:
            pass
        
        self.fastener_groups = []
        self.fastener_specs = {}
        self.other_elements = []
        
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)
        
        self.tab_widget = QTabWidget()
        self.main_layout.addWidget(self.tab_widget)
        
        self.calculation_tab = QWidget()
        self.groups_tab = QWidget()
        self.specs_tab = QWidget()
        self.help_tab = QWidget()
        
        self.tab_widget.addTab(self.calculation_tab, "Calculation")
        self.tab_widget.addTab(self.groups_tab, "Fastener Groups")
        self.tab_widget.addTab(self.specs_tab, "Fastener Specifications")
        self.tab_widget.addTab(self.help_tab, "Help and Documentation")
        
        self._setup_calculation_tab()
        self._setup_groups_tab()
        self._setup_specs_tab()
        self._setup_help_tab()
        
        self.worker = None
        
        self._load_fastener_specs()
        
        self.status_bar = self.statusBar()
        self.status_bar.setStyleSheet("""
            QStatusBar {
                background-color: #f0f0f0;
                border-top: 2px solid #cccccc;
                font-size: 10pt;
                color: #000000;
            }
            QStatusBar::item {
                border: none;
            }
        """)

        version = "v2.39"
        version_label = QLabel(version)
        version_label.setStyleSheet("color: #000000; font-size: 8pt;")
        self.status_bar.addWidget(version_label)

        author_label = QLabel("Ömer Kurt - Structural Analysis Engineer")
        author_label.setStyleSheet("color: #000000; font-size: 8pt;")
        self.status_bar.addPermanentWidget(author_label)

        self.show()
    
    def _setup_calculation_tab(self):
        layout = QVBoxLayout(self.calculation_tab)
        
        file_group = QGroupBox("File Selection")
        file_layout = QVBoxLayout()
        
        input_layout = QHBoxLayout()
        self.input_label = QLabel("Input BDF:")
        self.input_edit = QLineEdit()
        self.input_edit.setReadOnly(True)
        self.input_button = QPushButton("Browse...")
        self.input_button.clicked.connect(self._browse_input_file)
        input_layout.addWidget(self.input_label)
        input_layout.addWidget(self.input_edit)
        input_layout.addWidget(self.input_button)
        file_layout.addLayout(input_layout)
        
        sets_layout = QHBoxLayout()
        self.sets_label = QLabel("Fastener Sets BDF (optional):")
        self.sets_edit = QLineEdit()
        self.sets_edit.setReadOnly(True)
        self.sets_button = QPushButton("Browse...")
        self.sets_button.clicked.connect(self._browse_sets_file)
        self.parse_button = QPushButton("Parse BDF for Fastener Sets")
        self.parse_button.clicked.connect(self._parse_bdf)
        
        sets_layout.addWidget(self.sets_label)
        sets_layout.addWidget(self.sets_edit)
        sets_layout.addWidget(self.sets_button)
        sets_layout.addWidget(self.parse_button)
        file_layout.addLayout(sets_layout)
        
        output_layout = QHBoxLayout()
        self.output_label = QLabel("Output BDF:")
        self.output_edit = QLineEdit()
        self.output_edit.setReadOnly(True)
        self.output_button = QPushButton("Browse...")
        self.output_button.clicked.connect(self._browse_output_file)
        
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.output_edit)
        output_layout.addWidget(self.output_button)
        file_layout.addLayout(output_layout)
        
        sets_output_layout = QHBoxLayout()
        self.sets_output_label = QLabel("Fastener Sets Output BDF (optional):")
        self.sets_output_edit = QLineEdit()
        self.sets_output_edit.setReadOnly(True)
        self.sets_output_button = QPushButton("Browse...")
        self.sets_output_button.clicked.connect(self._browse_sets_output_file)
        self.export_sets_button = QPushButton("Export Updated Sets BDF")
        self.export_sets_button.clicked.connect(self._export_sets_bdf)
        
        sets_output_layout.addWidget(self.sets_output_label)
        sets_output_layout.addWidget(self.sets_output_edit)
        sets_output_layout.addWidget(self.sets_output_button)
        sets_output_layout.addWidget(self.export_sets_button)
        file_layout.addLayout(sets_output_layout)
        
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        # Create parallel processing section
        parallel_group = QGroupBox("Parallel Processing Settings")
        parallel_layout = QFormLayout()
        
        # Number of workers
        worker_layout = QHBoxLayout()
        self.num_workers_spin = QSpinBox()
        self.num_workers_spin.setRange(1, mp.cpu_count())
        self.num_workers_spin.setValue(min(4, mp.cpu_count()))
        self.num_workers_spin.setToolTip(f"Number of parallel workers (1-{mp.cpu_count()}, recommended: {min(4, mp.cpu_count())})")
        self.num_workers_spin.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        worker_layout.addWidget(self.num_workers_spin)

        # Worker info moved to right side
        worker_info = QLabel(f"Available CPU cores: {mp.cpu_count()}")
        worker_info.setStyleSheet("color: #666666; font-size: 9pt;")
        worker_layout.addWidget(worker_info)
        
        parallel_layout.addRow("Number of Workers:", worker_layout)
        
        parallel_group.setLayout(parallel_layout)
        layout.addWidget(parallel_group)
        
        # Connection Processing Mode
        connection_mode_group = QGroupBox("Connection Processing Mode")
        connection_mode_layout = QVBoxLayout()
        
        self.connection_mode_group = QButtonGroup()
        
        self.convert_to_airbus_radio = QRadioButton("Convert Multi-Plate Connections to Airbus Method")
        self.convert_to_airbus_radio.setToolTip(
            "Airbus Method:\n"
            "• Creates additional CBUSH spanning outer plates with only axial stiffness\n"
            "• Sets regular CBUSHes to have only shear stiffness (K1=100)\n"
            "• Only applies to multi-plate connections (3+ plates)\n"
            "• Converts existing normal connections to Airbus method"
        )
        
        self.convert_to_normal_radio = QRadioButton("Convert Airbus Connections to Normal Method")
        self.convert_to_normal_radio.setToolTip(
            "Normal Method:\n"
            "• Removes Airbus CBUSH elements and their properties\n"
            "• Updates regular CBUSHes with full stiffness calculations\n"
            "• Removes Airbus elements from fastener groups\n"
            "• Converts existing Airbus connections to normal method"
        )
        
        self.mix_mode_radio = QRadioButton("Preserve Existing Connection Types (Mix Mode Support)")
        self.mix_mode_radio.setToolTip(
            "Mix Mode:\n"
            "• Preserves existing Airbus connections as Airbus method\n"
            "• Preserves existing normal connections as normal method\n"
            "• Does not convert between connection types\n"
            "• Calculates stiffness appropriate for each connection type"
        )
        self.mix_mode_radio.setChecked(True)  # Default selection
        
        self.connection_mode_group.addButton(self.convert_to_airbus_radio, 0)
        self.connection_mode_group.addButton(self.convert_to_normal_radio, 1)
        self.connection_mode_group.addButton(self.mix_mode_radio, 2)
        
        connection_mode_layout.addWidget(self.convert_to_airbus_radio)
        connection_mode_layout.addWidget(self.convert_to_normal_radio)
        connection_mode_layout.addWidget(self.mix_mode_radio)
        
        connection_mode_group.setLayout(connection_mode_layout)
        layout.addWidget(connection_mode_group)
        
        # Create default fastener parameters section
        default_group = QGroupBox("Default Fastener Parameters")
        default_layout = QFormLayout()
        
        # Fastener specification with checkbox on the right
        spec_layout = QHBoxLayout()
        self.default_spec_combo = QComboBox()
        self.default_spec_combo.addItem("")
        self.default_spec_combo.currentTextChanged.connect(self._on_default_spec_changed)
        self.default_spec_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        spec_layout.addWidget(self.default_spec_combo)
        
        self.use_spec_check = QCheckBox("Use Fastener Specification")
        self.use_spec_check.stateChanged.connect(self._on_use_spec_changed)
        spec_layout.addWidget(self.use_spec_check)
        
        default_layout.addRow("Fastener Specification:", spec_layout)
        
        # Fastener diameter in mm (dropdown)
        self.diameter_label = QLabel("Fastener Diameter:")
        self.diameter_combo = QComboBox()
        self.diameter_combo.addItems(["4.78", "6.35", "7.92"])
        self.diameter_combo.setCurrentText("6.35")
        default_layout.addRow(self.diameter_label, self.diameter_combo)
        
        # Fastener type
        self.fastener_type_label = QLabel("Fastener Type:")
        self.fastener_type_combo = QComboBox()
        self.fastener_type_combo.addItems(["Bolt", "Rivet", "Hi-Lok", "Blind Rivet"])
        default_layout.addRow(self.fastener_type_label, self.fastener_type_combo)
        
        # Fastener material
        self.fastener_material_label = QLabel("Fastener Material:")
        self.fastener_material_combo = QComboBox()
        self.fastener_material_combo.addItems(["Titanium","Steel", "Aluminum", "Custom"])
        self.fastener_material_combo.currentTextChanged.connect(self._update_fastener_e)
        default_layout.addRow(self.fastener_material_label, self.fastener_material_combo)
        
        # Fastener elasticity in GPa
        self.fastener_e_label = QLabel("Elastic Modulus:")
        self.fastener_e_spin = QDoubleSpinBox()
        self.fastener_e_spin.setRange(10.0, 500.0)
        self.fastener_e_spin.setSingleStep(5.0)
        self.fastener_e_spin.setDecimals(1)
        self.fastener_e_spin.setValue(110.0)
        self.fastener_e_spin.setSuffix(" GPa")
        self.fastener_e_spin.setEnabled(False)
        default_layout.addRow(self.fastener_e_label, self.fastener_e_spin)
        
        # Material orientation approach
        self.material_orientation_label = QLabel("Material Orientation Approach:")
        self.material_orientation_combo = QComboBox()
        self.material_orientation_combo.addItems([
            "Directional Modulus (Ex, Ey from Jones Formula)",
            "Effective Modulus (E_eff = √(E1×E2))"
        ])
        self.material_orientation_combo.setToolTip(
            "Directional: Uses transformed Ex for K2, Ey for K3 based on material orientation\n"
            "Effective: Uses geometric mean of E1×E2 for both K2 and K3"
        )
        default_layout.addRow(self.material_orientation_label, self.material_orientation_combo)
        
        # K4, K5, K6 values
        self.k4_spin = QDoubleSpinBox()
        self.k4_spin.setRange(0.1, 1e10)
        self.k4_spin.setSingleStep(1.0)
        self.k4_spin.setDecimals(1)
        self.k4_spin.setValue(100.0)
        default_layout.addRow("K4 (Torsional Stiffness):", self.k4_spin)
        
        self.k5_spin = QDoubleSpinBox()
        self.k5_spin.setRange(1.0, 1e12)
        self.k5_spin.setSingleStep(1e7)
        self.k5_spin.setDecimals(0)
        self.k5_spin.setValue(1e8)
        default_layout.addRow("K5 (Rotational Stiffness 1):", self.k5_spin)
        
        self.k6_spin = QDoubleSpinBox()
        self.k6_spin.setRange(1.0, 1e12)
        self.k6_spin.setSingleStep(1e7)
        self.k6_spin.setDecimals(0)
        self.k6_spin.setValue(1e8)
        default_layout.addRow("K6 (Rotational Stiffness 2):", self.k6_spin)
        
        # Multi-plate connection method
        self.connection_method_label = QLabel("Multi-Plate Connection Method:")
        self.connection_method_combo = QComboBox()
        self.connection_method_combo.addItems([
            "Two-Plate Enforcement (Single Shear for All)",
            "HyperMesh Method (n=2 for Multi-Plate)",
            "Normal Multi-Plate Analysis"
        ])
        self.connection_method_combo.setToolTip(
            "Two-Plate: Forces all connections to be analyzed as two-plate (single shear)\n"
            "HyperMesh: Uses HyperMesh method - n=2 for multi-plate, outer/inner plate logic\n"
            "Normal: Uses actual multi-plate connection analysis"
        )
        default_layout.addRow(self.connection_method_label, self.connection_method_combo)
        
        default_group.setLayout(default_layout)
        layout.addWidget(default_group)
        
        # Other elements section
        other_group = QGroupBox("Other Elements to Update")
        other_layout = QVBoxLayout()
        
        other_info = QLabel("Define additional CBUSH element IDs not in fastener groups that should be updated with default parameters:")
        other_info.setWordWrap(True)
        other_layout.addWidget(other_info)
        
        self.other_elements_edit = QTextEdit()
        self.other_elements_edit.setMaximumHeight(35)
        self.other_elements_edit.setPlaceholderText("Enter comma-separated or space-separated element IDs (leave empty to update all if no fastener groups)")
        other_layout.addWidget(self.other_elements_edit)
        
        other_group.setLayout(other_layout)
        layout.addWidget(other_group)
        
        # Create a log area
        log_group = QGroupBox("Calculation Log")
        log_layout = QVBoxLayout()
        
        log_container = QHBoxLayout()
        
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(150)
        self.log_text.setText("Ready to calculate stiffness values using parallel processing...")
        log_container.addWidget(self.log_text)
        
        log_layout.addLayout(log_container)
        log_group.setLayout(log_layout)
        layout.addWidget(log_group)
        
        # Create progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)
        
        # Create calculate button
        self.calculate_button = QPushButton("Calculate and Update CBUSH Stiffness")
        self.calculate_button.setStyleSheet("font-weight: bold; padding: 10px;")
        self.calculate_button.clicked.connect(self._start_calculation)
        layout.addWidget(self.calculate_button)
        
        # Initialize UI state
        self._on_use_spec_changed()
    
    def _setup_groups_tab(self):
        """Setup the fastener groups tab UI"""
        layout = QVBoxLayout(self.groups_tab)
        
        # Instructions
        instructions = QLabel(
            "Fastener groups define which elements use specific fastener configurations. Groups can be "
            "automatically detected from a sets BDF file or created manually for specific elements."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)
        
        # Fastener groups table
        self.groups_table = QTableWidget(0, 6)
        self.groups_table.setHorizontalHeaderLabels(["Group Name", "Export Name", "Spec", "Diameter", "# Elements", "Issues"])
        self.groups_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.groups_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.groups_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.groups_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.groups_table.itemSelectionChanged.connect(self._update_group_details)
        
        layout.addWidget(self.groups_table)
        
        # Group details and element list
        details_splitter = QSplitter(Qt.Horizontal)
        
        # Group details
        details_widget = QWidget()
        details_layout = QVBoxLayout(details_widget)
        details_group = QGroupBox("Group Details")
        details_inner_layout = QVBoxLayout()
        self.group_details = QLabel("No group selected")
        self.group_details.setWordWrap(True)
        self.group_details.setAlignment(Qt.AlignTop)
        details_inner_layout.addWidget(self.group_details)
        details_group.setLayout(details_inner_layout)
        details_layout.addWidget(details_group)
        
        # Element list
        elements_widget = QWidget()
        elements_layout = QVBoxLayout(elements_widget)
        elements_group = QGroupBox("Element IDs")
        elements_inner_layout = QVBoxLayout()
        self.elements_list = QListWidget()
        elements_inner_layout.addWidget(self.elements_list)
        elements_group.setLayout(elements_inner_layout)
        elements_layout.addWidget(elements_group)
        
        details_splitter.addWidget(details_widget)
        details_splitter.addWidget(elements_widget)
        layout.addWidget(details_splitter)
        
        # Buttons for group management
        button_layout = QHBoxLayout()
        self.add_group_button = QPushButton("Add Manual Group")
        self.add_group_button.clicked.connect(self._add_group)
        self.edit_group_button = QPushButton("Edit Group")
        self.edit_group_button.clicked.connect(self._edit_group)
        self.delete_group_button = QPushButton("Delete Group")
        self.delete_group_button.clicked.connect(self._delete_group)
        
        button_layout.addWidget(self.add_group_button)
        button_layout.addWidget(self.edit_group_button)
        button_layout.addWidget(self.delete_group_button)
        
        layout.addLayout(button_layout)
    
    def _setup_specs_tab(self):
        """Setup the fastener specifications tab UI"""
        layout = QVBoxLayout(self.specs_tab)
        
        # Instructions
        instructions = QLabel(
            "Define fastener specifications for each fastener type detected in the BDF file. "
            "Each spec (like ELS438) needs to be mapped to a material and type for proper stiffness calculation."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)
        
        # Specifications table
        self.specs_table = QTableWidget(0, 4)
        self.specs_table.setHorizontalHeaderLabels(["Spec Name", "Fastener Type", "Material", "Elastic Modulus"])
        self.specs_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.specs_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.specs_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.specs_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        layout.addWidget(self.specs_table)
        
        # Buttons for spec management
        button_layout = QHBoxLayout()
        self.add_spec_button = QPushButton("Add Specification")
        self.add_spec_button.clicked.connect(self._add_specification)
        self.edit_spec_button = QPushButton("Edit Selected Specification")
        self.edit_spec_button.clicked.connect(self._edit_specification)
        self.delete_spec_button = QPushButton("Delete Selected Specification")
        self.delete_spec_button.clicked.connect(self._delete_specification)
        self.clear_specs_button = QPushButton("Clear All Specifications")
        self.clear_specs_button.clicked.connect(self._clear_specifications)
        
        button_layout.addWidget(self.add_spec_button)
        button_layout.addWidget(self.edit_spec_button)
        button_layout.addWidget(self.delete_spec_button)
        button_layout.addWidget(self.clear_specs_button)
        button_layout.addStretch()
        
        layout.addLayout(button_layout)
        
        # Import/Export buttons
        io_layout = QHBoxLayout()
        self.import_specs_button = QPushButton("Import Specifications")
        self.import_specs_button.clicked.connect(self._import_specifications)
        self.export_specs_button = QPushButton("Export Specifications")
        self.export_specs_button.clicked.connect(self._export_specifications)
        
        io_layout.addWidget(self.import_specs_button)
        io_layout.addWidget(self.export_specs_button)
        io_layout.addStretch()
        
        layout.addLayout(io_layout)
    
    def _setup_help_tab(self):
        """Setup the help tab UI with comprehensive documentation"""
        layout = QHBoxLayout(self.help_tab)
        
        # Create splitter for topic navigation and content
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)
        
        # Left side: Topic navigation
        nav_widget = QWidget()
        nav_layout = QVBoxLayout(nav_widget)
        nav_layout.setContentsMargins(5, 5, 5, 5)
        
        nav_label = QLabel("Documentation Topics")
        nav_label.setStyleSheet("font-weight: bold; font-size: 12pt; margin: 5px;")
        nav_layout.addWidget(nav_label)
        
        self.topic_list = QListWidget()
        self.topic_list.setMaximumWidth(250)
        self.topic_list.setStyleSheet("QListWidget { border: 1px solid #cccccc; }")
        
        # Add documentation topics
        topics = [
            ("Overview", "overview"),
            ("Huth Formula Theory", "huth_theory"),
            ("Material Orientation", "material_orientation"),
            ("Connection Types", "connection_types"),
            ("Airbus Method", "airbus_method"),
            ("Multi-Plate Analysis", "multiplate_analysis"),
            ("Fastener Specifications", "fastener_specs"),
            # ("Parallel Processing", "parallel_processing"),
            ("File Formats", "file_formats"),
            ("Coordinate Systems", "coordinate_systems"),
            # ("Classical Lamination Theory", "clt_theory"),
            # ("Validation & QA", "validation"),
            # ("Troubleshooting", "troubleshooting"),
            # ("Examples & Use Cases", "examples")
        ]
        
        for topic_name, topic_id in topics:
            item = QListWidgetItem(topic_name)
            item.setData(Qt.UserRole, topic_id)
            self.topic_list.addItem(item)
        
        self.topic_list.currentItemChanged.connect(self._on_topic_changed)
        nav_layout.addWidget(self.topic_list)
        
        # Right side: Content display
        self.help_content = QTextEdit()
        self.help_content.setReadOnly(True)
        self.help_content.setStyleSheet("QTextEdit { background-color: white; border: 1px solid #cccccc; }")
        
        splitter.addWidget(nav_widget)
        splitter.addWidget(self.help_content)
        splitter.setSizes([250, 800])
        
        # Set initial content
        self.topic_list.setCurrentRow(0)

    def _on_topic_changed(self, current, previous):
        """Handle topic selection change"""
        if current:
            topic_id = current.data(Qt.UserRole)
            content = self._get_topic_content(topic_id)
            self.help_content.setHtml(content)

    def _get_topic_content(self, topic_id):
        return get_topic_content(topic_id)
    def _on_use_spec_changed(self):
        """Handle use spec checkbox change"""
        use_spec = self.use_spec_check.isChecked()
        
        # Enable/disable spec combo
        self.default_spec_combo.setEnabled(use_spec)
        
        # Enable/disable individual parameter controls
        self.fastener_type_combo.setEnabled(not use_spec)
        self.fastener_material_combo.setEnabled(not use_spec)
        self.fastener_e_spin.setEnabled(not use_spec and self.fastener_material_combo.currentText() == "Custom")
        
        # Update specs combo
        self._update_default_spec_combo()
        
        if use_spec:
            self._on_default_spec_changed()
    
    def _on_default_spec_changed(self):
        """Handle default spec selection change"""
        if not self.use_spec_check.isChecked():
            return
            
        spec_name = self.default_spec_combo.currentText()
        if spec_name in self.fastener_specs:
            spec = self.fastener_specs[spec_name]
            # Update the UI to show spec values (read-only)
            self.fastener_type_combo.setCurrentText(spec.fastener_type)
            self.fastener_material_combo.setCurrentText(spec.material)
            self.fastener_e_spin.setValue(spec.elastic_modulus)
    
    def _update_default_spec_combo(self):
        """Update the default spec combo with available specs"""
        current = self.default_spec_combo.currentText()
        self.default_spec_combo.clear()
        self.default_spec_combo.addItem("")  # Empty option
        
        if self.fastener_specs:
            self.default_spec_combo.addItems(list(self.fastener_specs.keys()))
            if current in self.fastener_specs:
                self.default_spec_combo.setCurrentText(current)
    
    def _update_fastener_e(self, material):
        """Update fastener elastic modulus based on material selection"""
        if self.use_spec_check.isChecked():
            return  # Don't update if using spec
            
        if material == "Steel":
            self.fastener_e_spin.setValue(210.0)  # GPa
            self.fastener_e_spin.setEnabled(False)
        elif material == "Titanium":
            self.fastener_e_spin.setValue(110.0)  # GPa
            self.fastener_e_spin.setEnabled(False)
        elif material == "Aluminum":
            self.fastener_e_spin.setValue(70.0)   # GPa
            self.fastener_e_spin.setEnabled(False)
        else:  # Custom
            self.fastener_e_spin.setEnabled(True)
    
    def _browse_input_file(self):
        """Open file dialog to select input BDF file"""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select Input BDF File", "", "Nastran BDF Files (*.bdf *.dat *.nas);;All Files (*.*)")
        
        if filename:
            self.input_edit.setText(filename)
            
            # Auto-fill output filename
            if not self.output_edit.text():
                base, ext = os.path.splitext(filename)
                self.output_edit.setText(f"{base}_modified{ext}")
    
    def _browse_sets_file(self):
        """Open file dialog to select sets BDF file"""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select Fastener Sets BDF File", "", "Nastran BDF Files (*.bdf *.dat *.nas);;All Files (*.*)")
        
        if filename:
            self.sets_edit.setText(filename)
            
            # Auto-fill sets output filename
            if not self.sets_output_edit.text():
                base, ext = os.path.splitext(filename)
                self.sets_output_edit.setText(f"{base}_updated{ext}")
    
    def _browse_output_file(self):
        """Open file dialog to select output BDF file"""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Select Output BDF File", "", "Nastran BDF Files (*.bdf *.dat *.nas);;All Files (*.*)")
        
        if filename:
            self.output_edit.setText(filename)
    
    def _browse_sets_output_file(self):
        """Open file dialog to select sets output BDF file"""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Select Fastener Sets Output BDF File", "", "Nastran BDF Files (*.bdf *.dat *.nas);;All Files (*.*)")
        
        if filename:
            self.sets_output_edit.setText(filename)
    
    def _parse_bdf(self):
        """Parse BDF file to extract fastener sets"""
        sets_file = self.sets_edit.text()
        
        if not sets_file:
            QMessageBox.warning(self, "Input Error", "Please select a fastener sets BDF file.")
            return
        
        # Update log
        self._append_log("Parsing BDF file for fastener sets...")
        self.progress_bar.setValue(10)
        
        # Parse the BDF file
        parser = BdfParser(sets_file)
        success = parser.parse()
        
        if not success:
            QMessageBox.critical(self, "Parsing Error", "Failed to parse BDF file.")
            self._append_log("Error parsing BDF file.")
            self.progress_bar.setValue(0)
            return
        
        # Get fastener groups and update UI
        new_groups = parser.fastener_groups
        
        # Remove any existing groups from the same file (non-manual)
        self.fastener_groups = [g for g in self.fastener_groups if g.manual]
        
        # Add the new groups
        self.fastener_groups.extend(new_groups)
        self._update_fastener_groups_table()
        
        # Create default fastener specs if they don't exist
        for group in new_groups:
            if group.spec_name and group.spec_name not in self.fastener_specs:
                spec = FastenerSpec(group.spec_name)
                self.fastener_specs[group.spec_name] = spec
        
        # Update fastener specs table and combo
        self._update_fastener_specs_table()
        self._update_default_spec_combo()
        
        # Complete progress
        self.progress_bar.setValue(100)
        
        # Count different types of issues
        problematic_count = sum(1 for g in new_groups if g.problematic)
        diameter_missing_count = sum(1 for g in new_groups if g.diameter_missing)
        
        self._append_log(f"Parsed BDF file successfully. Found {len(new_groups)} fastener groups.")
        if problematic_count > 0:
            self._append_log(f"Warning: {problematic_count} groups had parsing issues.")
        if diameter_missing_count > 0:
            self._append_log(f"Warning: {diameter_missing_count} groups are missing diameter information.")
        
        # Switch to groups tab
        if new_groups:
            self.tab_widget.setCurrentWidget(self.groups_tab)
    
    def _export_sets_bdf(self):
        """Export fastener groups as BDF sets"""
        if not self.fastener_groups:
            QMessageBox.warning(self, "No Groups", "No fastener groups to export.")
            return
        
        # Use the sets output file if specified, otherwise prompt
        output_file = self.sets_output_edit.text()
        if not output_file:
            filename, _ = QFileDialog.getSaveFileName(
                self, "Export Sets BDF", "", "Nastran BDF Files (*.bdf *.dat *.nas);;All Files (*.*)")
            if not filename:
                return
            output_file = filename
        
        parser = BdfParser("")  # Empty filename since we're only exporting
        success = parser.export_sets_bdf(output_file, self.fastener_groups)
        
        if success:
            QMessageBox.information(self, "Export Successful", f"Sets exported to {output_file}")
            self._append_log(f"Exported {len(self.fastener_groups)} fastener groups to {output_file}")
        else:
            QMessageBox.critical(self, "Export Error", "Failed to export sets BDF file.")
    
    def _update_fastener_groups_table(self):
        """Update the fastener groups table"""
        # Clear existing rows
        self.groups_table.setRowCount(0)
        
        # Add rows for each fastener group
        self.groups_table.setRowCount(len(self.fastener_groups))
        
        for i, group in enumerate(self.fastener_groups):
            # Group name
            item = QTableWidgetItem(group.set_name)
            if group.problematic or group.diameter_missing:
                item.setBackground(QColor(255, 255, 200))  # Light yellow for issues
            self.groups_table.setItem(i, 0, item)
            
            # Export name (read-only)
            export_item = QTableWidgetItem(group.export_name)
            export_item.setFlags(export_item.flags() & ~Qt.ItemIsEditable)  # Make read-only
            export_item.setBackground(QColor(240, 240, 240))  # Light gray background
            self.groups_table.setItem(i, 1, export_item)
            
            # Fastener spec
            self.groups_table.setItem(i, 2, QTableWidgetItem(group.spec_name))
            
            # Diameter
            self.groups_table.setItem(i, 3, QTableWidgetItem(f"{group.diameter_mm:.2f} mm"))
            
            # Element count
            self.groups_table.setItem(i, 4, QTableWidgetItem(str(len(group.element_ids))))
            
            # Issues
            issues = []
            if group.diameter_missing:
                issues.append("Missing Diameter")
            if group.problematic:
                issues.append("Parsing Error")
            if not group.spec_name:
                issues.append("Missing Spec")
            if group.manual:
                issues.append("Manual")
            
            issues_text = ", ".join(issues) if issues else "None"
            self.groups_table.setItem(i, 5, QTableWidgetItem(issues_text))
    
    def _update_fastener_specs_table(self):
        """Update the fastener specifications table"""
        # Clear existing rows
        self.specs_table.setRowCount(0)
        
        # Add rows for each fastener spec
        self.specs_table.setRowCount(len(self.fastener_specs))
        
        i = 0
        for spec_name, spec in self.fastener_specs.items():
            # Spec name
            self.specs_table.setItem(i, 0, QTableWidgetItem(spec.spec_name))
            
            # Fastener type
            self.specs_table.setItem(i, 1, QTableWidgetItem(spec.fastener_type))
            
            # Material
            self.specs_table.setItem(i, 2, QTableWidgetItem(spec.material))
            
            # Elastic modulus
            self.specs_table.setItem(i, 3, QTableWidgetItem(f"{spec.elastic_modulus:.1f} GPa"))
            
            i += 1
    
    def _update_group_details(self):
        """Update the group details display"""
        current_row = self.groups_table.currentRow()
        if current_row >= 0 and current_row < len(self.fastener_groups):
            group = self.fastener_groups[current_row]
            
            # Create details text with proper line breaks
            details_lines = [
                f"Group Name: {group.set_name}",
                f"Export Name: {group.export_name}",
                f"Original Name: {group.original_name}",
                f"Fastener Specification: {group.spec_name or 'Not specified'}",
                f"Fastener Diameter: {group.diameter_mm:.2f} mm",
                f"Number of Elements: {len(group.element_ids)}",
                f"Source: {'Manual' if group.manual else 'SET BDF'}",
            ]
            
            # Add issue details
            issues = []
            if group.diameter_missing:
                issues.append("Missing diameter information")
            if group.problematic:
                issues.append("Parsing errors detected")
            if not group.spec_name:
                issues.append("No fastener specification")
            
            if issues:
                details_lines.append(f"Issues: {', '.join(issues)}")
            else:
                details_lines.append("Issues: None")
            
            self.group_details.setText('\n'.join(details_lines))
            
            # Update element list
            self.elements_list.clear()
            for elem_id in group.element_ids:
                self.elements_list.addItem(str(elem_id))
        else:
            self.group_details.setText("No group selected")
            self.elements_list.clear()
    
    def _add_group(self):
        """Add a new manual fastener group"""
        dialog = FastenerGroupDialog(specs=self.fastener_specs, parent=self)
        if dialog.exec_() == QDialog.Accepted:
            group = dialog.get_group()
            self.fastener_groups.append(group)
            self._update_fastener_groups_table()
    
    def _edit_group(self):
        """Edit the selected fastener group"""
        current_row = self.groups_table.currentRow()
        if current_row >= 0 and current_row < len(self.fastener_groups):
            group = self.fastener_groups[current_row]
                
            dialog = FastenerGroupDialog(group, self.fastener_specs, parent=self)
            if dialog.exec_() == QDialog.Accepted:
                updated_group = dialog.get_group()
                # Update the group in place to maintain reference
                self.fastener_groups[current_row] = updated_group
                self._update_fastener_groups_table()
                self._update_group_details()
    
    def _delete_group(self):
        """Delete the selected fastener group"""
        current_row = self.groups_table.currentRow()
        if current_row >= 0 and current_row < len(self.fastener_groups):
            if QMessageBox.question(self, "Confirm Deletion", 
                                    f"Delete group '{self.fastener_groups[current_row].set_name}'?",
                                    QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
                del self.fastener_groups[current_row]
                self._update_fastener_groups_table()
                self._update_group_details()
    
    def _add_specification(self):
        """Add a new fastener specification"""
        # Get a name from the user
        name, ok = QInputDialog.getText(self, "New Specification", "Enter specification name:")
        if ok and name:
            if name in self.fastener_specs:
                QMessageBox.warning(self, "Duplicate", f"Specification '{name}' already exists.")
                return
                
            spec = FastenerSpec(name)
            dialog = FastenerSpecDialog(spec, parent=self)
            if dialog.exec_() == QDialog.Accepted:
                self.fastener_specs[name] = dialog.get_spec()
                self._update_fastener_specs_table()
                self._update_default_spec_combo()
                self._save_fastener_specs()
                
                # Update any groups using this spec
                self._update_groups_with_spec(name)
    
    def _edit_specification(self):
        """Edit the selected fastener specification"""
        current_row = self.specs_table.currentRow()
        if current_row >= 0 and current_row < self.specs_table.rowCount():
            spec_name = self.specs_table.item(current_row, 0).text()
            
            if spec_name in self.fastener_specs:
                dialog = FastenerSpecDialog(self.fastener_specs[spec_name], parent=self)
                if dialog.exec_() == QDialog.Accepted:
                    self.fastener_specs[spec_name] = dialog.get_spec()
                    self._update_fastener_specs_table()
                    self._update_default_spec_combo()
                    self._save_fastener_specs()
                    
                    # Update any groups using this spec
                    self._update_groups_with_spec(spec_name)
    
    def _update_groups_with_spec(self, spec_name):
        """Update fastener groups that use the specified spec"""
        updated = False
        for group in self.fastener_groups:
            if group.spec_name == spec_name:
                group.update_export_name()  # Update export name after spec change
                if not group.manual:
                    group.update_set_name()
                updated = True
        
        if updated:
            self._update_fastener_groups_table()
            self._update_group_details()
    
    def _delete_specification(self):
        """Delete the selected fastener specification"""
        current_row = self.specs_table.currentRow()
        if current_row >= 0 and current_row < self.specs_table.rowCount():
            spec_name = self.specs_table.item(current_row, 0).text()
            
            # Check if spec is used by any fastener groups
            used_by = [g.set_name for g in self.fastener_groups if g.spec_name == spec_name]
            if used_by:
                QMessageBox.warning(self, "Cannot Delete", 
                                    f"Specification '{spec_name}' is used by fastener groups: {', '.join(used_by)}")
                return
                
            if QMessageBox.question(self, "Confirm Deletion", 
                                    f"Delete specification '{spec_name}'?",
                                    QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
                del self.fastener_specs[spec_name]
                self._update_fastener_specs_table()
                self._update_default_spec_combo()
                self._save_fastener_specs()
    
    def _clear_specifications(self):
        """Clear all fastener specifications"""
        # Check if any specs are used
        used_specs = set(g.spec_name for g in self.fastener_groups if g.spec_name)
        if used_specs:
            QMessageBox.warning(self, "Cannot Clear", 
                                f"Some specifications are in use: {', '.join(used_specs)}")
            return
            
        if QMessageBox.question(self, "Confirm Clear", 
                                "Clear all fastener specifications?",
                                QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
            self.fastener_specs.clear()
            self._update_fastener_specs_table()
            self._update_default_spec_combo()
            self._save_fastener_specs()
    
    def _import_specifications(self):
        """Import fastener specifications from file"""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Import Fastener Specifications", "", "JSON Files (*.json);;All Files (*.*)")
        
        if filename:
            try:
                with open(filename, 'r') as f:
                    data = json.load(f)
                
                imported_count = 0
                for spec_name, spec_data in data.items():
                    if spec_name not in self.fastener_specs:
                        spec = FastenerSpec.from_dict(spec_data)
                        self.fastener_specs[spec_name] = spec
                        imported_count += 1
                
                self._update_fastener_specs_table()
                self._update_default_spec_combo()
                self._save_fastener_specs()
                
                QMessageBox.information(self, "Import Successful", 
                                        f"Imported {imported_count} fastener specifications.")
                self._append_log(f"Imported {imported_count} fastener specifications from {filename}")
                
            except Exception as e:
                QMessageBox.critical(self, "Import Error", f"Failed to import specifications: {str(e)}")
    
    def _export_specifications(self):
        """Export fastener specifications to file"""
        if not self.fastener_specs:
            QMessageBox.warning(self, "No Specifications", "No fastener specifications to export.")
            return
            
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Fastener Specifications", "fastener_specs.json", "JSON Files (*.json);;All Files (*.*)")
        
        if filename:
            try:
                data = {name: spec.to_dict() for name, spec in self.fastener_specs.items()}
                with open(filename, 'w') as f:
                    json.dump(data, f, indent=2)
                
                QMessageBox.information(self, "Export Successful", 
                                        f"Exported {len(self.fastener_specs)} fastener specifications.")
                self._append_log(f"Exported {len(self.fastener_specs)} fastener specifications to {filename}")
                
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export specifications: {str(e)}")
    
    def _load_fastener_specs(self):
        """Load fastener specs from default file if it exists"""
        # Get the directory where the script/exe is located
        if getattr(sys, 'frozen', False):
            # Running as exe
            app_dir = os.path.dirname(sys.executable)
        else:
            # Running as script
            app_dir = os.path.dirname(os.path.abspath(__file__))
        
        default_file = os.path.join(app_dir, "fastener_specs.json")
        
        if os.path.exists(default_file):
            try:
                with open(default_file, 'r') as f:
                    data = json.load(f)
                
                for spec_name, spec_data in data.items():
                    spec = FastenerSpec.from_dict(spec_data)
                    self.fastener_specs[spec_name] = spec
                
                self._update_fastener_specs_table()
                self._update_default_spec_combo()
                
            except Exception as e:
                print(f"Error loading default fastener specs: {str(e)}")
    
    def _save_fastener_specs(self):
        """Save fastener specs to default file"""
        # Get the directory where the script/exe is located
        if getattr(sys, 'frozen', False):
            # Running as exe
            app_dir = os.path.dirname(sys.executable)
        else:
            # Running as script
            app_dir = os.path.dirname(os.path.abspath(__file__))
        
        default_file = os.path.join(app_dir, "fastener_specs.json")
        
        try:
            data = {name: spec.to_dict() for name, spec in self.fastener_specs.items()}
            with open(default_file, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            print(f"Error saving fastener specs: {str(e)}")
    
    def _get_connection_mode(self):
        """Get the selected connection mode"""
        checked_button = self.connection_mode_group.checkedButton()
        if checked_button == self.convert_to_airbus_radio:
            return 'convert_to_airbus'
        elif checked_button == self.convert_to_normal_radio:
            return 'convert_to_normal'
        else:  # mix_mode_radio
            return 'mix_mode'
    
    def _start_calculation(self):
        """Start the parallel calculation process"""
        # Validate inputs
        input_file = self.input_edit.text()
        output_file = self.output_edit.text()
        
        if not input_file:
            QMessageBox.warning(self, "Input Error", "Please select an input BDF file.")
            return
        
        if not output_file:
            QMessageBox.warning(self, "Input Error", "Please select an output BDF file.")
            return
        
        # Parse other elements
        other_elements = []
        other_text = self.other_elements_edit.toPlainText().strip()
        if other_text:
            for item in other_text.replace(',', ' ').split():
                try:
                    other_elements.append(int(item))
                except ValueError:
                    pass
        
        # Get number of workers
        num_workers = self.num_workers_spin.value()
        
        # Get connection mode
        connection_mode = self._get_connection_mode()
        
        # Disable UI during calculation
        self._set_ui_enabled(False)
        
        # Reset progress bar and log
        self.progress_bar.setValue(0)
        self.log_text.clear()
        self._append_log(f"Starting parallel calculation with {num_workers} workers...")
        self._append_log(f"Connection mode: {connection_mode}")
        
        # Get connection method
        if "Two-Plate" in self.connection_method_combo.currentText():
            connection_method = "two_plate"
        elif "HyperMesh" in self.connection_method_combo.currentText():
            connection_method = "hypermesh"
        else:
            connection_method = "normal"
        
        # Get material orientation approach
        use_material_orientation = "Directional" in self.material_orientation_combo.currentText()
        
        # Get default parameters
        if self.use_spec_check.isChecked() and self.default_spec_combo.currentText():
            # Use selected fastener spec
            spec_name = self.default_spec_combo.currentText()
            if spec_name in self.fastener_specs:
                spec = self.fastener_specs[spec_name]
                default_parameters = {
                    'fastener_diameter_mm': float(self.diameter_combo.currentText()),
                    'fastener_type': spec.fastener_type.lower(),
                    'fastener_e': spec.elastic_modulus,  # GPa
                    'connection_method': connection_method,
                    'k4': self.k4_spin.value(),
                    'k5': self.k5_spin.value(),
                    'k6': self.k6_spin.value(),
                    'use_material_orientation': use_material_orientation
                }
            else:
                QMessageBox.warning(self, "Invalid Specification", "Selected fastener specification not found.")
                self._set_ui_enabled(True)
                return
        else:
            # Use individual parameters
            default_parameters = {
                'fastener_diameter_mm': float(self.diameter_combo.currentText()),
                'fastener_type': self.fastener_type_combo.currentText().lower(),
                'fastener_e': self.fastener_e_spin.value(),  # GPa
                'connection_method': connection_method,
                'k4': self.k4_spin.value(),
                'k5': self.k5_spin.value(),
                'k6': self.k6_spin.value(),
                'use_material_orientation': use_material_orientation
            }
        
        # Create and start parallel worker thread
        self.worker = ParallelCalculationWorker(
            input_file, output_file, default_parameters, self.fastener_groups, 
            self.fastener_specs, other_elements, num_workers, connection_mode)
        self.worker.progress_updated.connect(self.progress_bar.setValue)
        self.worker.calculation_finished.connect(self._on_calculation_finished)
        self.worker.log_message.connect(self._append_log)
        self.worker.groups_updated.connect(self._on_groups_updated)  # New signal connection
        self.worker.start()
    
    def _on_groups_updated(self):
        """Handle groups updated signal from worker"""
        # Update the GUI tables to reflect changes
        self._update_fastener_groups_table()
        self._update_group_details()
        self._append_log("GUI updated to reflect group changes")
    
    def _append_log(self, message):
        """Append message to log"""
        self.log_text.append(message)
        # Auto-scroll to bottom
        cursor = self.log_text.textCursor()
        cursor.movePosition(cursor.End)
        self.log_text.setTextCursor(cursor)
    
    def _on_calculation_finished(self, success, message):
        """Handle when calculation is finished"""
        # Re-enable UI
        self._set_ui_enabled(True)
        
        # Update fastener groups table after calculation completion
        self._update_fastener_groups_table()
        self._update_group_details()
        
        # Show message
        if success:
            QMessageBox.information(self, "Parallel Calculation Complete", message)
        else:
            QMessageBox.critical(self, "Calculation Error", message)
    
    def _set_ui_enabled(self, enabled):
        """Enable or disable UI elements during calculation"""
        # File selection
        self.input_button.setEnabled(enabled)
        self.sets_button.setEnabled(enabled)
        self.output_button.setEnabled(enabled)
        self.sets_output_button.setEnabled(enabled)
        self.parse_button.setEnabled(enabled)
        self.export_sets_button.setEnabled(enabled)
        self.calculate_button.setEnabled(enabled)
        
        # Parallel processing settings
        self.num_workers_spin.setEnabled(enabled)
        
        # Calculation tab
        self.use_spec_check.setEnabled(enabled)
        self.default_spec_combo.setEnabled(enabled and self.use_spec_check.isChecked())
        self.diameter_combo.setEnabled(enabled)
        self.connection_method_combo.setEnabled(enabled)
        self.material_orientation_combo.setEnabled(enabled)
        
        # Connection mode radio buttons
        self.convert_to_airbus_radio.setEnabled(enabled)
        self.convert_to_normal_radio.setEnabled(enabled)
        self.mix_mode_radio.setEnabled(enabled)
        
        self.fastener_type_combo.setEnabled(enabled and not self.use_spec_check.isChecked())
        self.fastener_material_combo.setEnabled(enabled and not self.use_spec_check.isChecked())
        self.fastener_e_spin.setEnabled(enabled and not self.use_spec_check.isChecked() and 
                                       self.fastener_material_combo.currentText() == "Custom")
        self.k4_spin.setEnabled(enabled)
        self.k5_spin.setEnabled(enabled)
        self.k6_spin.setEnabled(enabled)
        self.other_elements_edit.setEnabled(enabled)
        
        # Groups tab
        self.add_group_button.setEnabled(enabled)
        self.edit_group_button.setEnabled(enabled)
        self.delete_group_button.setEnabled(enabled)
        self.groups_table.setEnabled(enabled)
        
        # Specs tab
        self.add_spec_button.setEnabled(enabled)
        self.edit_spec_button.setEnabled(enabled)
        self.delete_spec_button.setEnabled(enabled)
        self.clear_specs_button.setEnabled(enabled)
        self.import_specs_button.setEnabled(enabled)
        self.export_specs_button.setEnabled(enabled)
        self.specs_table.setEnabled(enabled)


def main():
    """Main function to start the application"""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("CBUSH Huth Stiffness Updater")
    app.setApplicationVersion("2.39")
    
    window = HuthCalculatorGUI()
    sys.exit(app.exec_())


if __name__ == "__main__":
    mp.freeze_support()
    main()
