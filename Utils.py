import numpy as np
from pyNastran.bdf.bdf import BDF
import os
import tempfile

def truncate_status(text, limit=85):
    """Truncate status text if it exceeds the limit, adding '...' at the end."""
    if len(text) <= limit:
        return text
    return text[:limit-3] + "..."
def parse_column_format(column_name):
    """
    Parse column header to extract name and format information.
    Format should be specified after last hyphen with underscore separator.
    Examples: 
    - "Thickness-of-Webs-Float_2" -> ("Thickness of Webs", float, 2)
    - "Material-Int" -> ("Material", int, None)
    - "Type-Str" -> ("Type", str, None)
    If no format specified, infers from data type.
    """
    parts = column_name.split('-')
    if len(parts)>1:
        format_part = parts[-1].lower()
        name = ' '.join(parts[:-1])
    else:
        format_part= []
        name=column_name
    # Check if format is specified in last part
        
    if 'float' in format_part:
        try:
            decimals = int(format_part.split('_')[1])
            return name, float, decimals
        except (IndexError, ValueError):
            return name, float, 2  # Default 2 decimals
    elif 'int' in format_part:
        return name, int, None
    elif 'str' in format_part:
        try:
            decimals = int(format_part.split('_')[1])
            return name, str, decimals
        except (IndexError, ValueError):
            return name, str, 3  # Default 3 decimals
    else:
        # No format specified, return original name
        return column_name, None, None

def prepare_bdf_with_bulk_section(bdf_path):
    """
    Prepare BDF file by ensuring it has proper BULK section structure.
    
    Args:
        bdf_path (str): Path to the original BDF file
        
    Returns:
        str: Path to the prepared temporary BDF file
    """
    try:
        with open(bdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except UnicodeDecodeError:
        # Try with different encoding if UTF-8 fails
        with open(bdf_path, 'r', encoding='latin-1', errors='ignore') as f:
            content = f.read()
    
    lines = content.splitlines()
    
    # Find BEGIN BULK and ENDDATA indices
    begin_bulk_idx = None
    enddata_idx = None
    
    for i, line in enumerate(lines):
        line_upper = line.strip().upper()
        if line_upper.startswith('BEGIN BULK'):
            begin_bulk_idx = i
        elif line_upper.startswith('ENDDATA'):
            enddata_idx = i
            break  # ENDDATA should be the last, so break here
    
    # Extract the bulk data content
    if begin_bulk_idx is not None and enddata_idx is not None:
        # Get content between BEGIN BULK and ENDDATA
        bulk_content = lines[begin_bulk_idx + 1:enddata_idx]
    elif begin_bulk_idx is not None:
        # Has BEGIN BULK but no ENDDATA - take everything after BEGIN BULK
        bulk_content = lines[begin_bulk_idx + 1:]
    elif enddata_idx is not None:
        # Has ENDDATA but no BEGIN BULK - take everything before ENDDATA
        bulk_content = lines[:enddata_idx]
    else:
        # No BEGIN BULK or ENDDATA - assume entire file is bulk data
        bulk_content = lines
    
    # Filter out comment lines (starting with $) and empty lines
    filtered_bulk_content = []
    for line in bulk_content:
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('$'):
            filtered_bulk_content.append(line)
    
    # Create the properly formatted BDF content
    formatted_content = ['CEND', 'BEGIN BULK'] + filtered_bulk_content + ['ENDDATA']
    
    # Create temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix='.bdf', prefix='temp_bdf_')
    try:
        with os.fdopen(temp_fd, 'w') as temp_file:
            temp_file.write('\n'.join(formatted_content))
        return temp_path
    except:
        os.close(temp_fd)
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise

def extract_bulk_data_only(bdf_path):
    """
    Extract only the bulk data content (between BEGIN BULK and ENDDATA) without headers.
    
    Args:
        bdf_path (str): Path to the original BDF file
        
    Returns:
        str: Path to the temporary BDF file with only bulk data
    """
    try:
        with open(bdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(bdf_path, 'r', encoding='latin-1', errors='ignore') as f:
            content = f.read()
    
    lines = content.splitlines()
    
    # Find BEGIN BULK and ENDDATA indices
    begin_bulk_idx = None
    enddata_idx = None
    
    for i, line in enumerate(lines):
        line_upper = line.strip().upper()
        if line_upper.startswith('BEGIN BULK'):
            begin_bulk_idx = i
        elif line_upper.startswith('ENDDATA'):
            enddata_idx = i
            break
    
    # Extract only the bulk data content
    if begin_bulk_idx is not None and enddata_idx is not None:
        bulk_content = lines[begin_bulk_idx + 1:enddata_idx]
    elif begin_bulk_idx is not None:
        bulk_content = lines[begin_bulk_idx + 1:]
    elif enddata_idx is not None:
        bulk_content = lines[:enddata_idx]
    else:
        bulk_content = lines

    # Filter out comment lines (starting with $) and empty lines
    filtered_bulk_content = []
    for line in bulk_content:
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('$'):
            filtered_bulk_content.append(line)
    # Create temporary file with only bulk data
    temp_fd, temp_path = tempfile.mkstemp(suffix='.bdf', prefix='bulk_only_')
    try:
        with os.fdopen(temp_fd, 'w') as temp_file:
            temp_file.write('\n'.join(filtered_bulk_content))
        return temp_path
    except:
        os.close(temp_fd)
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise
def extract_bulk_data_only_grid(bdf_path):
    """
    Extract only the bulk data content (between BEGIN BULK and ENDDATA) without headers.
    
    Args:
        bdf_path (str): Path to the original BDF file
        
    Returns:
        str: Path to the temporary BDF file with only bulk data
    """
    try:
        with open(bdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(bdf_path, 'r', encoding='latin-1', errors='ignore') as f:
            content = f.read()
    
    lines = content.splitlines()
    
    # Find BEGIN BULK and ENDDATA indices
    begin_bulk_idx = None
    enddata_idx = None
    
    for i, line in enumerate(lines):
        line_upper = line.strip().upper()
        if line_upper.startswith('GRID'):
            begin_bulk_idx = i
        elif line_upper.startswith('ENDDATA'):
            enddata_idx = i
            break
    
    # Extract only the bulk data content
    if begin_bulk_idx is not None and enddata_idx is not None:
        bulk_content = lines[begin_bulk_idx:enddata_idx]
    elif begin_bulk_idx is not None:
        bulk_content = lines[begin_bulk_idx:]
    elif enddata_idx is not None:
        bulk_content = lines[:enddata_idx]
    else:
        bulk_content = lines

    # Filter out comment lines (starting with $) and empty lines
    filtered_bulk_content = []
    for line in bulk_content:
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('$'):
            filtered_bulk_content.append(line)
    # Create temporary file with only bulk data
    temp_fd, temp_path = tempfile.mkstemp(suffix='.bdf', prefix='grid_bulk_only_')
    try:
        with os.fdopen(temp_fd, 'w') as temp_file:
            temp_file.write('\n'.join(filtered_bulk_content))
        return temp_path
    except:
        os.close(temp_fd)
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise
def try_bdf_file(bdf_path):
    """
    Robust BDF file reading function that tries multiple approaches to successfully read the file.
    
    Args:
        bdf_path (str): Path to the BDF file
        
    Returns:
        BDF: Successfully loaded pyNastran BDF model
        
    Raises:
        Exception: If all reading attempts fail
    """
    model = BDF()
    temp_files_to_cleanup = []
    
    try:
        # Method 1: Try reading directly with different combinations
        read_methods = [
            {'punch': False, 'xref': True},
            {'punch': True, 'xref': True},
            {'punch': False, 'xref': False},
            {'punch': True, 'xref': False}
        ]
        
        for method in read_methods:
            try:
                print(f"Trying to read BDF with punch={method['punch']}, xref={method['xref']}")
                model.read_bdf(bdf_path, **method)
                print("Successfully read BDF file directly")
                return model
            except Exception as e:
                print(f"Direct read failed with punch={method['punch']}, xref={method['xref']}: {str(e)}")
                continue
        
        # Method 2: Try with properly formatted BULK section (xref=True, punch=False)
        try:
            print("Trying to read with properly formatted BULK section...")
            formatted_path = prepare_bdf_with_bulk_section(bdf_path)
            temp_files_to_cleanup.append(formatted_path)
            
            model.read_bdf(formatted_path, punch=False, xref=True, validate=False)
            print("Successfully read BDF with formatted BULK section")
            return model
            
        except Exception as e:
            print(f"Formatted BULK section read failed: {str(e)}")
        
        # Method 3: Try with bulk data only (xref=False, punch=True)
        try:
            print("Trying to read with bulk data only...")
            bulk_only_path = extract_bulk_data_only(bdf_path)
            temp_files_to_cleanup.append(bulk_only_path)
            
            model.read_bdf(bulk_only_path, punch=True, xref=False,validate=False)
            print("Successfully read BDF with bulk data only")
            return model
            
        except Exception as e:
            print(f"Bulk data only read failed: {str(e)}")
        
        # Method 4: Try with bulk data only start from grid (xref=False, punch=True)
        try:
            print("Trying to read with bulk data only (from GRID start)...")
            bulk_only_grid_path = extract_bulk_data_only_grid(bdf_path)
            temp_files_to_cleanup.append(bulk_only_grid_path)
            
            model.read_bdf(bulk_only_grid_path, punch=True, xref=False,validate=False)
            print("Successfully read BDF with bulk data only from GRID start")
            return model
                
        except Exception as e:
            print(f"Bulk data only read failed: {str(e)}")
        # If all methods fail, raise the last exception
        raise Exception("All BDF reading methods failed. The file may be corrupted or have missing required data.")
        
    finally:
        # Clean up temporary files
        for temp_file in temp_files_to_cleanup:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            except:
                pass  # Ignore cleanup errors

def find_free_edges(shell_elements, cbar_elements, nodes,offset_show):
    """
    Find free/boundary edges for both shell elements and expanded CBAR elements.
    
    Args:
        shell_elements (dict): Dictionary of shell elements
        cbar_elements (dict): Dictionary of CBAR elements
        nodes (dict): Original node coordinates
        adjusted_nodes (dict): Adjusted node coordinates (after transformation)
        
    Returns:
        list: List of edge coordinate pairs for all boundary edges
    """
    coords = np.array(list(nodes.values()))
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    
    # Adjust node coordinates
    adjusted_nodes = {nid: coords - min_coords for nid, coords in nodes.items()}
    
    # Group elements by property
    property_groups = {}
    for eid, elem_data in shell_elements.items():
        pid = elem_data['pid']
        if pid not in property_groups:
            property_groups[pid] = []
        property_groups[pid].append((eid, elem_data['nodes']))

    # Find edges that only appear once in each property group
    shell_free_edges = []
    for pid, elems in property_groups.items():
        edge_count = {}
        for eid, element_nodes in elems:
            for i in range(len(element_nodes)):
                edge = tuple(sorted((element_nodes[i], 
                                   element_nodes[(i + 1) % len(element_nodes)])))
                edge_count[edge] = edge_count.get(edge, 0) + 1

        # Only keep edges that appear once (boundary edges)
        free_edges = [edge for edge, count in edge_count.items() if count == 1]
        shell_free_edges.extend(free_edges)

    # Create edge coordinates for shell elements using adjusted coordinates
    edge_coords = []
    for edge in shell_free_edges:
        coords = [adjusted_nodes[node_id] for node_id in edge]
        edge_coords.append(coords)

    # Handle CBAR elements
    for eid, elem_data in cbar_elements.items():
        # Get width and orientation info
        width = elem_data['width']
        orientation = np.array(elem_data['orientation'])
        orientation = orientation / np.linalg.norm(orientation)
        direction = orientation * (width/2)

        # Get adjusted node positions
        start_node = adjusted_nodes[elem_data['nodes'][0]]
        end_node = adjusted_nodes[elem_data['nodes'][1]]
        
        # Apply offsets if they exist
        if 'offsets' in elem_data and offset_show:
            offsets = elem_data['offsets']
            offt = elem_data['offt']
            
            # Apply start offset if it exists
            if 'wa' in offsets:
                # Apply offset based on OFFT flag
                if offt[1] == 'G':  # Global coordinate system
                    start_node = start_node + offsets['wa']
                elif offt[1] == 'B':  # Basic coordinate system
                    # Transform offset to local
                    # This would need coordinate system transformation logic pass for now
                    start_node = start_node + offsets['wa']  # Simplified
            
            # Apply end offset if it exists
            if 'wb' in offsets:
                if offt[2] == 'G':  # Global coordinate system
                    end_node = end_node + offsets['wb']
                elif offt[2] == 'B':  # Basic coordinate system
                    # This would need coordinate system transformation logic pass for now
                    end_node = end_node + offsets['wb']  # Simplified
                    
        # Calculate the four corners of the expanded CBAR
        v1 = start_node + direction  # top front
        v2 = start_node - direction  # bottom front
        v3 = end_node - direction    # bottom back
        v4 = end_node + direction    # top back

        # Add the side edges (perpendicular to orientation)
        edge_coords.extend([
            [v1, v4],  # front edge
            [v2, v3]   # back edge
        ])

    # Convert edge coordinates to the format expected by plotly
    plotly_edges = []
    for edge in edge_coords:
        plotly_edges.extend([edge[0], edge[1], [None, None, None]])

    return np.array(plotly_edges)
def extract_original_mesh_edges(shell_elements, cbar_elements, nodes, offset_show=False):
    """
    Extract original mesh edges from shell and CBAR elements.
    
    Args:
        shell_elements (dict): Shell element data from BDF
        cbar_elements (dict): CBAR element data from BDF
        nodes (dict): Node coordinates
        offset_show (bool): Whether to apply offsets to CBAR elements
        
    Returns:
        np.array: Array of edge coordinates for plotting
    """
    coords = np.array(list(nodes.values()))
    min_coords = np.min(coords, axis=0)
    
    # # Adjust node coordinates
    adjusted_nodes = {nid: coords - min_coords for nid, coords in nodes.items()}
    # Initialize list to store all edge segments
    edge_coords = []
    
    # Process shell elements (quads and triangles)
    for eid, elem_data in shell_elements.items():
        elem_nodes = elem_data['nodes']
        
        # For each edge in the element
        for i in range(len(elem_nodes)):
            node1 = elem_nodes[i]
            node2 = elem_nodes[(i + 1) % len(elem_nodes)]
            
            # Add edge coordinates
            if node1 in adjusted_nodes and node2 in adjusted_nodes:
                edge_coords.append(adjusted_nodes[node1])
                edge_coords.append(adjusted_nodes[node2])
                edge_coords.append([None, None, None])  # Separator for plotly
    
    # Process CBAR elements
    for eid, elem_data in cbar_elements.items():
        # CBAR elements have 2 nodes (start and end)
        node1, node2 = elem_data['nodes']
        
        if node1 in adjusted_nodes and node2 in adjusted_nodes:
            # Apply offsets if they exist
            start_node = adjusted_nodes[node1].copy()
            end_node = adjusted_nodes[node2].copy()
            
            if offset_show and 'offsets' in elem_data:
                offsets = elem_data['offsets']
                offt = elem_data['offt']
                
                # Apply start offset if it exists
                if 'wa' in offsets:
                    if offt[1] == 'G':  # Global coordinate system
                        start_node = start_node + offsets['wa']
                    elif offt[1] == 'B':  # Basic coordinate system
                        start_node = start_node + offsets['wa']  # Simplified
                
                # Apply end offset if it exists
                if 'wb' in offsets:
                    if offt[2] == 'G':  # Global coordinate system
                        end_node = end_node + offsets['wb']
                    elif offt[2] == 'B':  # Basic coordinate system
                        end_node = end_node + offsets['wb']  # Simplified
            
            # Add edge coordinates
            edge_coords.append(start_node)
            edge_coords.append(end_node)
            edge_coords.append([None, None, None])  # Separator for plotly
    
    return np.array(edge_coords)