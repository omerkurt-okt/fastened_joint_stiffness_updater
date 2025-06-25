
def get_topic_content(topic_id):
    """Get content for specified topic"""
    
    base_style = """
    <style>
    body { font-family: 'Segoe UI', Arial, sans-serif; line-height: 1.6; color: #333; margin: 20px; }
    h1 { color: #2c5aa0; border-bottom: 2px solid #2c5aa0; padding-bottom: 5px; }
    h2 { color: #2c5aa0; margin-top: 25px; margin-bottom: 10px; }
    h3 { color: #4a90e2; margin-top: 20px; margin-bottom: 8px; }
    .formula { background-color: #f8f9fa; border-left: 4px solid #2c5aa0; padding: 10px; margin: 10px 0; font-family: 'Courier New', monospace; }
    .note { background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; margin: 10px 0; border-radius: 4px; }
    .warning { background-color: #f8d7da; border: 1px solid #f5c6cb; padding: 10px; margin: 10px 0; border-radius: 4px; }
    .code { background-color: #f8f9fa; padding: 2px 4px; font-family: 'Courier New', monospace; border-radius: 3px; }
    table { border-collapse: collapse; width: 100%; margin: 10px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #f2f2f2; font-weight: bold; }
    ul, ol { margin: 10px 0; padding-left: 30px; }
    </style>
    """
    
    if topic_id == "overview":
        return base_style + """
        <h1>CBUSH Huth Stiffness Updater - Overview</h1>
        
        <h2>Purpose</h2>
        <p>The CBUSH Huth Stiffness Updater is an advanced finite element analysis utility designed to calculate and update CBUSH element stiffness values using the Huth formula for fastened joints.</p>
        
        <h2>Key Capabilities</h2>
        <ul>
            <li><strong>Huth Formula Implementation:</strong> Accurate calculation of fastener stiffness based on plate properties and fastener characteristics</li>
            <li><strong>Material Orientation Support:</strong> Advanced handling of composite materials with proper coordinate transformations for Ex and Ey calculations</li>
            <li><strong>Multi-Plate Connection Analysis:</strong> Automatic identification and analysis of complex multi-plate joints</li>
            <li><strong>Airbus Method Support:</strong> Implementation of Airbus-specific fastener modeling approach</li>
            <li><strong>Parallel Processing:</strong> Multi-core processing for large models with thousands of fasteners</li>
            <li><strong>Classical Lamination Theory:</strong> Full CLT implementation for composite plates (PCOMP, PCOMPG, PCOMPP) property calculations</li>
            <li><strong>Fastener Database:</strong> Comprehensive fastener specification management</li>
            <li><strong>Set Management:</strong> Import/export of fastener groups from Nastran SET cards</li>
        </ul>
        
        <h2>Workflow</h2>
        <ol>
            <li><strong>Model Preparation:</strong> Load Nastran BDF file containing CBUSH elements</li>
            <li><strong>Fastener Configuration:</strong> Define fastener groups and specifications</li>
            <li><strong>Connection Analysis:</strong> Automatic identification of multi-plate connections</li>
            <li><strong>Stiffness Calculation:</strong> Apply Huth formula with material orientation</li>
            <li><strong>Model Update:</strong> Update PBUSH properties with calculated stiffness values</li>
            <li><strong>Output Generation:</strong> Export modified BDF and updated fastener sets</li>
        </ol>

        """
    
    elif topic_id == "huth_theory":
        return base_style + """
        <h1>Huth Formula Theory</h1>
        
        <h2>Background</h2>
        <p>The Huth formula is an empirical method for calculating the shear stiffness of fastened joints, widely used in aerospace structural analysis. It accounts for the flexibility of both the fastener and the connected plates.</p>
        
        <h2>Fundamental Equations</h2>
        
        <h3>Shear Stiffness Calculation</h3>
        <div class="formula">
        K_shear = 1 / C_total<br>
        <br>
        Where:<br>
        C_total = [(t₁ + t₂)/(2d)]^a × [b₁×C₁ + b₂×C₂]<br>
        
        C₁ = 1/(t₁×E₁) + 1/(2×t₁×E_f)<br>
        C₂ = 1/(t₂×E₂) + 1/(2×t₂×E_f)
        </div>
        
        <h3>Variable Definitions</h3>
        <table>
        <tr><th>Variable</th><th>Definition</th><th>Units</th></tr>
        <tr><td>t₁, t₂</td><td>Thickness of plates 1 and 2</td><td>mm</td></tr>
        <tr><td>d</td><td>Fastener diameter</td><td>mm</td></tr>
        <tr><td>E₁, E₂</td><td>Elastic modulus of plates 1 and 2</td><td>MPa</td></tr>
        <tr><td>E_f</td><td>Fastener elastic modulus</td><td>MPa</td></tr>
        <tr><td>a</td><td>Geometric parameter (material dependent)</td><td>-</td></tr>
        <tr><td>b₁, b₂</td><td>Material parameters for plates 1 and 2</td><td>-</td></tr>
        </table>
        
        <h3>Material Parameters</h3>
        <table>
        <tr><th>Material Type</th><th>a</th><th>b</th></tr>
        <tr><td>Metal (Bolt)</td><td>2/3</td><td>3.0</td></tr>
        <tr><td>Metal (Rivet)</td><td>2/5</td><td>2.2</td></tr>
        <tr><td>Composite</td><td>2/3</td><td>4.2</td></tr>
        </table>
        
        <h3>Axial Stiffness Calculation</h3>
        <div class="formula">
        K_axial = (E_f × π × d²) / (4 × t_total)<br>
        <br>
        Where:<br>
        t_total = sum of all plate thicknesses in the connection
        </div>
        
        <h2>Multi-Plate Extensions</h2>
        
        <h3>Normal Multi-Plate Method</h3>
        <p>For connections with more than 2 plates, the method identifies the strongest middle plate and treats other plates as equivalent outer plates.</p>
        
        <h3>HyperMesh Method</h3>
        <p>Uses n=2 for multi-plate connections and applies specific logic for outer/inner plate identification:</p>
        <div class="formula">
        b₁ = b₁_base / n<br>
        b₂ = b₂_base / n²<br>
        <br>
        Where n = 2 for multi-plate connections
        </div>
        
        <h2>Directional Stiffness</h2>
        <p>When material orientation is considered, separate stiffness values are calculated for each fastener direction:</p>
        
        <div class="formula">
        K₂ = f(E_x_transformed)  // Fastener e₂ direction<br>
        K₃ = f(E_y_transformed)  // Fastener e₃ direction
        </div>
        
        <div class="note">
        <strong>Implementation Note:</strong> The tool automatically selects appropriate material parameters based on fastener type and plate materials. Composite materials use enhanced parameters to account for their anisotropic behavior.
        </div>
        """
    
    elif topic_id == "material_orientation":
        return base_style + """
        <h1>Material Orientation and Coordinate Transformations</h1>
        
        <h2>Overview</h2>
        <p>Material orientation significantly affects fastener stiffness in composite and orthotropic materials. This tool implements comprehensive coordinate transformation to account for material direction relative to fastener orientation.</p>
        
        <h2>Coordinate Systems Involved</h2>
        
        <h3>1. Global Coordinate System</h3>
        <p>The overall model coordinate system where all geometry is defined.</p>
        
        <h3>2. Element Coordinate System</h3>
        <p>Local coordinate system for each plate element:</p>
        <ul>
            <li><strong>CTRIA3:</strong> X-axis from G1 to G2 (side 1-2)</li>
            <li><strong>CQUAD4:</strong> Simplified to G1-G2 direction</li>
            <li><strong>Z-axis:</strong> Normal to plate surface (right-hand rule)</li>
        </ul>
        
        <h3>3. Material Coordinate System</h3>
        <p>Principal material directions (fiber directions for composites):</p>
        <ul>
            <li>Defined by <span class="code">theta_mcid</span> parameter on element cards</li>
            <li>Can be angle (degrees) or coordinate system ID</li>
            <li>Rotated relative to element coordinate system</li>
        </ul>
        
        <h3>4. Fastener Coordinate System</h3>
        <p>Local coordinate system for CBUSH element:</p>
        <ul>
            <li><strong>e₁:</strong> Along fastener axis (axial direction)</li>
            <li><strong>e₂:</strong> Fastener Y-direction (shear, K₂)</li>
            <li><strong>e₃:</strong> Fastener Z-direction (shear, K₃)</li>
        </ul>
        
        <h2>Transformation Process</h2>
        
        <h3>Step 1: Element Geometry Analysis</h3>
        <div class="formula">
        element_x = normalize(G2_position - G1_position)<br>
        plate_normal = normalize(cross(v1, v2))
        </div>
        
        <h3>Step 2: Material Direction Calculation</h3>
        <p>Depending on <span class="code">theta_mcid</span> type:</p>
        
        <h4>Direct Angle (degrees)</h4>
        <div class="formula">
        material_dir = rotate_in_plane(element_x, plate_normal, theta_rad)
        </div>
        
        <h4>Coordinate System Reference</h4>
        <div class="formula">
        cs_x_direction = coordinate_system[MCID].x_axis<br>
        material_dir = project_onto_plane(cs_x_direction, plate_normal)
        </div>
        
        <h3>Step 3: Fastener Coordinate System</h3>
        <div class="formula">
        fastener_e1 = normalize(node2_pos - node1_pos)<br>
        fastener_e3 = normalize(cross(fastener_e1, orientation_vector))<br>
        fastener_e2 = cross(fastener_e3, fastener_e1)
        </div>
        
        <h3>Step 4: Material-Fastener Angle Calculation</h3>
        <div class="formula">
        fastener_e2_on_plate = project_onto_plane(fastener_e2, plate_normal)<br>
        fastener_e3_on_plate = project_onto_plane(fastener_e3, plate_normal)<br>
        <br>
        theta_e2 = angle_between(material_dir, fastener_e2_on_plate)<br>
        theta_e3 = angle_between(material_dir, fastener_e3_on_plate)<br>
        </div>
        
        <h2>Jones Formula Transformation</h2>
        
        <h3>For Orthotropic Materials (MAT8)</h3>
        <div class="formula">
        1/E(θ) = (1/E₁)cos⁴θ + (1/E₂)sin⁴θ + (1/G₁₂ - 2ν₁₂/E₁)cos²θsin²θ<br>
        <br>
        Where:<br>
        θ = angle between material direction and loading direction<br>
        E₁, E₂ = principal elastic moduli<br>
        G₁₂ = in-plane shear modulus<br>
        ν₁₂ = major Poisson's ratio<br>
        </div>
        
        <h3>Application to Fastener Directions</h3>
        <div class="formula">
        E_e2 = jones_transform(E₁, E₂, G₁₂, ν₁₂, theta_e2)<br>
        E_e3 = jones_transform(E₁, E₂, G₁₂, ν₁₂, theta_e3)<br>
        <br>
        K₂ = huth_formula(E_e2, other_params)<br>
        K₃ = huth_formula(E_e3, other_params)<br>
        </div>
        
        <h2>Material Orientation Approaches</h2>
        
        <h3>Effective Modulus Approach</h3>
        <ul>
            <li>Uses geometric mean: E_eff = √(E₁ × E₂)</li>
            <li>Same stiffness for both K₂ and K₃</li>
            <li>Conservative and simple</li>
            <li>Suitable when orientation is not critical</li>
        </ul>
        
        <h3>Directional Modulus Approach</h3>
        <ul>
            <li>Calculates separate E_e2 and E_e3 using Jones formula</li>
            <li>Different stiffness for K₂ and K₃</li>
            <li>More accurate for oriented materials</li>
            <li>Requires proper element orientation definition</li>
        </ul>
        
        <div class="warning">
        <strong>Important:</strong> Material orientation transformation is applied only during stiffness calculation, not during property extraction. This prevents double transformation and ensures accuracy.
        </div>
        """
    
    elif topic_id == "connection_types":
        return base_style + """
        <h1>Connection Types and Analysis Methods</h1>
        
        <h2>Connection Identification</h2>
        <p>The tool automatically identifies multi-plate connections by finding CBUSH elements that share nodes. Connected CBUSHes are grouped into single connections for comprehensive analysis.</p>
        
        <h3>Connection Detection Algorithm</h3>
        <div class="formula">
        1. Start with unprocessed CBUSH element<br>
        2. Find all CBUSHes sharing nodes with current element<br>
        3. Recursively find CBUSHes connected to found elements<br>
        4. Group all connected CBUSHes into single connection<br>
        5. Repeat until all CBUSHes are processed
        </div>
        
        <h2>Connection Types</h2>
        
        <h3>Two-Plate Connections</h3>
        <p><strong>Characteristics:</strong></p>
        <ul>
            <li>Single CBUSH element</li>
            <li>Connects exactly two plates</li>
            <li>Standard single-shear joint</li>
        </ul>
        <p><strong>Analysis:</strong> Direct application of Huth formula with plate properties.</p>
        
        <h3>Multi-Plate Connections</h3>
        <p><strong>Characteristics:</strong></p>
        <ul>
            <li>Multiple CBUSH elements sharing nodes</li>
            <li>Three or more plates in stack-up</li>
            <li>Complex load transfer through fastener</li>
        </ul>
        <p><strong>Analysis Methods:</strong></p>
        
        <h4>1. Normal Multi-Plate Analysis</h4>
        <ul>
            <li>Identifies strongest middle plate</li>
            <li>Creates equivalent outer plate properties</li>
            <li>Applies modified Huth formula with n=2</li>
        </ul>
        
        <h4>2. HyperMesh Method</h4>
        <ul>
            <li>Uses specific outer/inner plate logic</li>
            <li>Applies n=2 scaling for multi-plate effects</li>
            <li>Compatible with HyperMesh fastener modeling</li>
        </ul>
        
        <h4>3. Two-Plate Enforcement</h4>
        <ul>
            <li>Forces all connections to be analyzed as two-plate</li>
            <li>Uses n=1 (single shear) for all connections</li>
            <li>Conservative approach</li>
        </ul>
        
        <h2>Plate Ordering and Properties</h2>
        
        <h3>Automatic Plate Ordering</h3>
        <p>For multi-plate connections, plates are automatically ordered from one end to the other:</p>
        <div class="formula">
        1. Build node-to-CBUSH mapping<br>
        2. Find end nodes (connected to only one CBUSH)<br>
        3. Traverse from one end to other following CBUSH chain<br>
        4. Order plates based on traversal sequence
        </div>
        
        <h3>Plate Property Extraction</h3>
        <p>Comprehensive support for various plate property types:</p>
        
        <table>
        <tr><th>Property Type</th><th>Material Types</th><th>Features</th></tr>
        <tr><td>PSHELL</td><td>MAT1, MAT8</td><td>Simple plate properties</td></tr>
        <tr><td>PCOMP</td><td>MAT8 layers</td><td>Classical lamination theory</td></tr>
        <tr><td>PCOMPG</td><td>Global plies</td><td>Global ply database</td></tr>
        <tr><td>PCOMPP</td><td>Parametric</td><td>Parametric composite definition</td></tr>
        </table>
        
        <h2>RBE3 Connection Handling</h2>
        <p>The tool properly handles CBUSH-to-plate connections through RBE3 elements:</p>
        
        <h3>Connection Path</h3>
        <div class="formula">
        CBUSH node → RBE3 dependent node → RBE3 independent nodes → Plate elements
        </div>
        
        <h3>Plate Identification</h3>
        <ol>
            <li>Check if CBUSH node is RBE3 dependent node</li>
            <li>Find RBE3 independent nodes</li>
            <li>Locate plate elements containing independent nodes</li>
            <li>Extract plate properties and thickness</li>
        </ol>
        
        <h2>Connection Method Selection Guide</h2>
        
        <table>
        <tr><th>Method</th><th>Best For</th><th>Characteristics</th></tr>
        <tr><td>Normal Multi-Plate</td><td>General analysis</td><td>Accurate multi-plate modeling</td></tr>
        <tr><td>HyperMesh Method</td><td>HyperMesh compatibility</td><td>Industry standard approach</td></tr>
        <tr><td>Two-Plate Enforcement</td><td>Conservative analysis</td><td>Simple, conservative results</td></tr>
        </table>
        
        <div class="note">
        <strong>Recommendation:</strong> Use Normal Multi-Plate for most applications. Use HyperMesh Method when compatibility with HyperMesh results is required. Use Two-Plate Enforcement for preliminary or conservative analysis.
        </div>
        """
    
    elif topic_id == "airbus_method":
        return base_style + """
        <h1>Airbus Method for Multi-Plate Connections</h1>
        
        <h2>Concept</h2>
        <p>The Airbus method is an alternative approach for modeling multi-plate fastened connections where the fastener's load transfer is split between axial and shear components using separate CBUSH elements.</p>
        
        <h2>Method Overview</h2>
        
        <h3>Dual CBUSH Approach</h3>
        <ul>
            <li><strong>Airbus CBUSH:</strong> Spans outer plates, carries only axial loads (K₁ >> K₂, K₃)</li>
            <li><strong>Regular CBUSHes:</strong> Between adjacent plates, carry only shear loads (K₁ ≈ 0)</li>
        </ul>
        
        <h3>Load Transfer Mechanism</h3>
        <div class="formula">
        Total Load = Axial Component + Shear Components<br>
        <br>
        Axial: Transferred through Airbus CBUSH (outer-to-outer)<br>
        Shear: Transferred through regular CBUSHes (plate-to-plate)
        </div>
        
        <h2>Stiffness Distribution</h2>
        
        <h3>Airbus CBUSH Properties</h3>
        <div class="formula">
        K₁ = (E_f × π × d²) / (4 × t_total)    // Full axial stiffness<br>
        K₂ = 100.0                            // Minimal shear stiffness<br>
        K₃ = 100.0                            // Minimal shear stiffness<br>
        K₄ = user_defined<br>
        K₅ = user_defined<br>
        K₆ = user_defined
        </div>
        
        <h3>Regular CBUSH Properties</h3>
        <div class="formula">
        K₁ = 100.0                            // Minimal axial stiffness<br>
        K₂ = Huth_formula(plate_props)        // Full shear stiffness<br>
        K₃ = Huth_formula(plate_props)        // Full shear stiffness<br>
        K₄ = user_defined<br>
        K₅ = user_defined<br>
        K₆ = user_defined
        </div>
        
        <h2>Airbus CBUSH Detection</h2>
        <p>Existing Airbus CBUSHes are automatically detected based on stiffness patterns:</p>
        
        <h3>Detection Criteria</h3>
        <div class="formula">
        IF (K₁ > 1000) AND (K₂ < 1000) AND (K₃ < 1000) AND (K₁ > 10 × max(K₂, K₃))<br>
        THEN element is classified as Airbus CBUSH
        </div>
        
        <h2>Connection Processing Modes</h2>
        
        <h3>Convert to Airbus Method</h3>
        <ul>
            <li>Creates new Airbus CBUSH spanning outer plates</li>
            <li>Modifies regular CBUSHes to carry only shear</li>
            <li>Adds new elements to appropriate fastener groups</li>
        </ul>
        
        <h3>Convert to Normal Method</h3>
        <ul>
            <li>Removes existing Airbus CBUSH elements</li>
            <li>Updates regular CBUSHes with full stiffness calculation</li>
            <li>Removes Airbus elements from fastener groups</li>
        </ul>
        
        <h3>Mix Mode (Preserve Existing)</h3>
        <ul>
            <li>Maintains existing connection types</li>
            <li>Calculates appropriate stiffness for each type</li>
            <li>No conversion between methods</li>
        </ul>
        
        <h2>Airbus CBUSH Creation Process</h2>
        
        <h3>Step 1: Identify Outer Nodes</h3>
        <div class="formula">
        1. Order plates from one end to other<br>
        2. Select outermost plates (first and last)<br>
        3. Find CBUSH nodes connected to outer plates<br>
        4. Select appropriate nodes for Airbus CBUSH
        </div>
        
        <h3>Step 2: Create New Elements</h3>
        <div class="formula">
        new_cbush_id = max(existing_element_ids) + 1<br>
        new_pbush_id = max(existing_property_ids) + 1<br>
        <br>
        CBUSH: new_cbush_id, new_pbush_id, [outer_node1, outer_node2]<br>
        PBUSH: new_pbush_id, [K₁_axial, 100, 100, K₄, K₅, K₆]
        </div>
        
        <h3>Step 3: Update Connection Data</h3>
        <ul>
            <li>Add new CBUSH to connection's Airbus element list</li>
            <li>Update connection classification</li>
            <li>Assign to appropriate fastener groups</li>
        </ul>
        
        <h2>Advantages and Limitations</h2>
        
        <h3>Advantages</h3>
        <ul>
            <li>Clear separation of axial and shear load paths</li>
            <li>Better representation of actual fastener behavior</li>
            <li>Improved numerical stability for some analyses</li>
            <li>Industry acceptance in certain applications</li>
        </ul>
        
        <h3>Limitations</h3>
        <ul>
            <li>Increased model complexity (more elements)</li>
            <li>Requires careful validation of results</li>
            <li>May not be suitable for all analysis types</li>
        </ul>
        
        <h2>When to Use Airbus Method</h2>
        
        <table>
        <tr><th>Use Case</th><th>Recommendation</th><th>Reason</th></tr>
        <tr><td>Multi-plate (3+ plates)</td><td>Consider Airbus</td><td>Better load path representation</td></tr>
        <tr><td>Dynamic analysis</td><td>Evaluate carefully</td><td>May affect modal characteristics</td></tr>
        </table>
        
        <div class="warning">
        <strong>Important:</strong> Always validate Airbus method results against validated/tested results. The method should be used consistently throughout the model for similar connections.
        </div>
        """
    
    elif topic_id == "multiplate_analysis":
        return base_style + """
        <h1>Multi-Plate Connection Analysis</h1>
        
        <h2>Introduction</h2>
        <p>Multi-plate connections involve three or more plates connected by a single fastener. These connections require special treatment to account for the complex load transfer and varying stiffness contributions of different plates.</p>
        
        <h2>Multi-Plate Challenges</h2>
        
        <h3>Load Transfer Complexity</h3>
        <ul>
            <li>Non-uniform load distribution between plates</li>
            <li>Plate sequence affects stiffness</li>
            <li>Material property variations</li>
            <li>Thickness variations</li>
        </ul>
        
        <h3>Modeling Considerations</h3>
        <ul>
            <li>Which plates to use in Huth formula</li>
            <li>How to handle middle plates</li>
            <li>Effective thickness calculation</li>
            <li>Material property averaging
            </ul>
       
       <h2>Analysis Methods Implemented</h2>
       
       <h3>1. Normal Multi-Plate Method</h3>
       <p>This method identifies the strongest middle plate and creates equivalent outer plate properties.</p>
       
       <h4>Algorithm</h4>
       <div class="formula">
       1. Order plates from one end to other<br>
       2. Calculate strength measure for each plate:<br>
       >   Strength = thickness × √(E/1000)<br>
       3. Identify strongest plate as "middle plate"<br>
       4. Average properties of remaining plates as "outer equivalent"<br>
       5. Apply Huth formula with n=2
       </div>
       
       <h4>Strength Calculation</h4>
       <div class="formula">
       Plate_Strength = t_plate × √(E_plate / 1000)<br>
       <br>
       Where:<br>
       t_plate = plate thickness (mm)<br>
       E_plate = plate elastic modulus (MPa)<br>
       </div>
       
       <h3>2. HyperMesh Method</h3>
       <p>Compatible with HyperMesh fastener modeling approach. Considers two plate at a time, goes CBUSH by CBUSH. Defines outer and middle plate between these two plates according to the following algorithm.</p>
       
       <h4>Plate Selection Logic</h4>
       <div class="formula">
       EVALUATE EACH CBUSH
       IF (any plate is outer plate):<br>
       >    middle_plate = non-outer plate<br>
       >    outer_plate = outer plate<br>
       ELSE:<br>
       >    middle_plate = stronger of the two plates<br>
       >    outer_plate = weaker of the two plates<br>
       </div>
       
       <h4>Multi-Plate Scaling</h4>
       <div class="formula">
       b₁ = b₁_base / n<br>
       b₂ = b₂_base / n²<br>
       <br>
       Where n = 2 for multi-plate connections
       </div>
       
       <h3>3. Two-Plate Enforcement</h3>
       <p>Simplifies all connections to two-plate analysis regardless of actual plate count.</p>
       
       <h4>Implementation</h4>
       <ul>
           <li>Uses first two plates found in connection</li>
           <li>Applies standard Huth formula with n=1</li>
           <li>Conservative approach</li>
           <li>Simplest to validate</li>
       </ul>
       
       <h2>Plate Property Averaging</h2>
       
       <h3>Thickness-Weighted Averaging</h3>
       <div class="formula">
       E_avg = Σ(E_i × t_i) / Σ(t_i)<br>
       t_avg = Σ(t_i) / n_plates<br>
       <br>
       Where:<br>
       E_i = elastic modulus of plate i<br>
       t_i = thickness of plate i<br>
       n_plates = number of plates being averaged
       </div>
       
       <h3>Geometric Mean (Alternative)</h3>
       <div class="formula">
       E_avg = (Π E_i)^(1/n)<br>
       <br>
       Used for highly dissimilar materials
       </div>
       
       <h2>Connection Strength Assessment</h2>
       
       <h3>Relative Strength Calculation</h3>
       <p>The tool calculates relative plate strength to determine which plate should be considered the "middle" plate in multi-plate connections:</p>
       
       <div class="formula">
       Strength_Factor = thickness × √(E_modulus / reference_modulus)<br>
       <br>
       Where reference_modulus = 1000 MPa (normalization factor)
       </div>
       
       <h3>Material Type Considerations</h3>
       <table>
       <tr><th>Material Type</th><th>Strength Calculation</th><th>Notes</th></tr>
       <tr><td>Isotropic (MAT1)</td><td>t × √(E/1000)</td><td>Direct calculation</td></tr>
       <tr><td>Orthotropic (MAT8)</td><td>t × √(E_effective OR E_directional/1000)</td><td>Uses geometric mean or directional modulus</td></tr>
       </table>
       
       <h2>Example Multi-Plate Analysis</h2>
       
       <h3>Three-Plate Stack-up</h3>
       <p>Consider a connection with three plates:</p>
       <ul>
           <li>Plate 1 (outer): Aluminum, 2.0mm, E=70000 MPa</li>
           <li>Plate 2 (middle): Steel, 1.5mm, E=210000 MPa</li>
           <li>Plate 3 (outer): Aluminum, 2.0mm, E=70000 MPa</li>
       </ul>
       
       <h4>Strength Calculation</h4>
       <div class="formula">
       Strength₁ = 2.0 × √(70000/1000) = 2.0 × 8.37 = 16.74<br>
       Strength₂ = 1.5 × √(210000/1000) = 1.5 × 14.49 = 21.74<br>
       Strength₃ = 2.0 × √(70000/1000) = 2.0 × 8.37 = 16.74<br>
       <br>
       → Plate 2 (steel) is identified as strongest (middle plate)
       </div>
       
       <h4>Outer Plate Equivalent</h4>
       <div class="formula">
       E_outer_equiv = (70000×2.0 + 70000×2.0) / (2.0 + 2.0) = 70000 MPa<br>
       t_outer_equiv = (2.0 + 2.0) / 2 = 2.0 mm<br>
       </div>
       
       <h4>Huth Formula Application</h4>
       <div class="formula">
       Middle plate: Steel, 1.5mm, 210000 MPa<br>
       Outer equivalent: Aluminum, 2.0mm, 70000 MPa<br>
       n = 2 (multi-plate factor)
       </div>
       
       <h2>Special Cases</h2>
       
       <h3>Identical Plates</h3>
       <p>When all plates have identical properties:</p>
       <ul>
           <li>No strength-based selection needed/li>
           <li>Apply appropriate multi-plate factors</li>
       </ul>
       
       <h3>Composite-Metal Combinations</h3>
       <p>Mixed material connections require careful treatment:</p>
       <ul>
           <li>Different material parameters (b values)</li>
           <li>Anisotropic vs isotropic behavior</li>
       </ul>
       
       <h3>Large Thickness Variations</h3>
       <p>When plates have very different thicknesses:</p>
       <ul>
           <li>Thickness-weighted averaging becomes critical</li>
           <li>Consider using Airbus method</li>
           <li>Validate with detailed analysis</li>
       </ul>
       
       <h2>Validation Recommendations</h2>
       
       <h3>Model Validation</h3>
       <ol>
           <li>Validate load paths and deformation patterns</li>
           <li>Compare with detailed 3D solid models if available</li>
           <li>Check against test data if available</li>
           <li>Perform sensitivity studies on key parameters</li>
       </ol>
       
       <h3>Results Checking</h3>
       <ul>
           <li>Verify stiffness values are physically reasonable</li>
           <li>Check for consistency across similar connections</li>
       </ul>
       """
   
    elif topic_id == "fastener_specs":
       return base_style + """
       <h1>Fastener Specifications and Database Management</h1>
       
       <h2>Overview</h2>
       <p>The fastener specification system allows users to define and manage different fastener types with their material properties and characteristics. This enables consistent application of fastener properties across large models.</p>
       
       <h2>Fastener Specification Components</h2>
       
       <h3>Basic Properties</h3>
       <table>
       <tr><th>Property</th><th>Description</th><th>Options/Range</th></tr>
       <tr><td>Specification Name</td><td>Unique identifier (e.g., "ELS438")</td><td>User-defined string</td></tr>
       <tr><td>Fastener Type</td><td>Mechanical fastener type</td><td>Bolt, Rivet, Hi-Lok, Blind Rivet</td></tr>
       <tr><td>Material</td><td>Fastener material</td><td>Titanium, Steel, Aluminum, Custom</td></tr>
       <tr><td>Elastic Modulus</td><td>Material elastic modulus</td><td>10-500 GPa</td></tr>
       </table>
       
       <h3>Predefined Material Properties</h3>
       <table>
       <tr><th>Material</th><th>Elastic Modulus (GPa)</th><th>Typical Applications</th></tr>
       <tr><td>Titanium</td><td>110</td><td>High-performance</td></tr>
       <tr><td>Steel</td><td>210</td><td>High strength</td></tr>
       <tr><td>Aluminum</td><td>70</td><td>Lightweight applications</td></tr>
       <tr><td>Custom</td><td>User-defined</td><td>Special materials</td></tr>
       </table>
       
       <h2>Fastener Type Effects on Huth Parameters</h2>
       
       <h3>Geometric Parameter (a)</h3>
       <table>
       <tr><th>Fastener Type</th><th>Parameter 'a'</th><th>Physical Meaning</th></tr>
       <tr><td>Bolt</td><td>2/3</td><td>Threaded fastener with high contact pressure</td></tr>
       <tr><td>Rivet</td><td>2/5</td><td>Deformed fastener with different bearing characteristics</td></tr>
       <tr><td>Hi-Lok</td><td>2/3</td><td>Similar to bolt characteristics</td></tr>
       <tr><td>Blind Rivet</td><td>2/5</td><td>Similar to rivet characteristics</td></tr>
       </table>
       
       <h3>Material Parameters (b)</h3>
       <p>The 'b' parameter in Huth formula is affected by both fastener type and plate material:</p>
       
       <div class="formula">
       b_effective = b_base_fastener × material_factor<br>
       <br>
       Where:<br>
       b_base_fastener = 3.0 (bolt/hi-lok) or 2.2 (rivet/blind rivet)<br>
       material_factor = 1.0 (metal) or 1.4 (composite)
       </div>
       
       <h2>Specification Management</h2>
       
       <h3>Creating Specifications</h3>
       <ol>
           <li>Navigate to "Fastener Specifications" tab</li>
           <li>Click "Add Specification"</li>
           <li>Enter unique specification name</li>
           <li>Select fastener type and material</li>
           <li>Adjust elastic modulus if needed</li>
           <li>Save specification</li>
       </ol>
       
       <h3>Applying Specifications</h3>
       <p>Specifications can be applied in two ways:</p>
       
       <h4>1. Default Parameters</h4>
       <ul>
           <li>Set default specification for all unassigned elements</li>
           <li>Use "Use Fastener Specification" checkbox</li>
           <li>Select specification from dropdown</li>
       </ul>
       
       <h4>2. Fastener Groups</h4>
       <ul>
           <li>Assign different specifications to different groups</li>
           <li>Automatic detection from SET names</li>
           <li>Manual assignment through group editing</li>
       </ul>
       
       <h2>SET-Based Specification Detection</h2>
       
       <h3>Naming Convention</h3>
       <p>The tool automatically detects fastener specifications from SET names using this pattern:</p>
       
       <div class="formula">
       xxx__xxx__DIAMETER__SPECIFICATION__xxx___xxx___xxx
       
       Example: "Frame200__Duct__6.35__ELS438__2__2__0"
       </div>
       
       <h3>Parsing Rules</h3>
       <ul>
           <li>Parts separated by double underscore "__"</li>
           <li>Part 1: General description, connected part 1</li>
           <li>Part 2: General description, connected part 2</li>
           <li>Part 3: Diameter in mm</li>
           <li>Part 4: Specification name</li>
           <li>Part 5: Row number (optional)</li>
           <li>Part 6: Connected plate number (optional)</li>
           <li>Part 7: Enumerator (optional)</li>
       </ul>
       
       <h3>Error Handling</h3>
       <p>When SET parsing encounters issues:</p>
       <ul>
           <li>Missing diameter: Group marked with "diameter_missing" flag</li>
           <li>Parsing errors: Group marked with "problematic" flag</li>
           <li>Unknown specifications: Default properties used with warning</li>
       </ul>
       
       <h2>Database Import/Export</h2>
       
       <h3>Export Format (JSON)</h3>
       <div class="formula">
       {<br>
         "ELS438": {<br>
           "spec_name": "ELS438",<br>
           "fastener_type": "Bolt",<br>
           "material": "Titanium",<br>
           "elastic_modulus": 110.0<br>
         },<br>
         "MS20470": {<br>
           "spec_name": "MS20470",<br>
           "fastener_type": "Rivet",<br>
           "material": "Aluminum",<br>
           "elastic_modulus": 70.0<br>
         }<br>
       }<br>
       </div>
       
       <h3>Persistent Storage</h3>
       <p>Specifications are automatically saved to:</p>
       <ul>
           <li>File: <span class="code">fastener_specs.json</span></li>
           <li>Location: Same directory as application</li>
           <li>Format: JSON with human-readable structure</li>
           <li>Auto-loaded: On application startup</li>
       </ul>
       
       
       <h2>Quality Assurance</h2>
       
       <h3>Specification Validation</h3>
       <ul>
           <li>Unique names enforced</li>
           <li>Reasonable elastic modulus ranges</li>
           <li>Consistent material-property relationships</li>
           <li>Warning for unusual combinations</li>
       </ul>
       
       <h3>Usage Tracking</h3>
       <p>The tool tracks specification usage:</p>
       <ul>
           <li>Shows which groups use each specification</li>
           <li>Prevents deletion of specifications in use</li>
           <li>Updates group names when specifications change</li>
           <li>Maintains consistency across model</li>
       </ul>
       
       <div class="note">
       <strong>Best Practice:</strong> Establish a consistent specification database for your team and share fastener_specs.json with your team-mates to use the same database.
       </div>
       """
   
    elif topic_id == "file_formats":
       return base_style + """
       <h1>File Formats and Data Exchange</h1>
       
       <h2>Supported File Formats</h2>
       
       <h3>Input Files</h3>
       <table>
       <tr><th>File Type</th><th>Extensions</th><th>Purpose</th><th>Required</th></tr>
       <tr><td>Nastran BDF</td><td>.bdf</td><td>Main finite element model</td><td>Yes</td></tr>
       <tr><td>Fastener Sets BDF</td><td>.bdf</td><td>SET definitions for fastener groups with HM comments</td><td>Optional</td></tr>
       <tr><td>Fastener Specs JSON</td><td>.json</td><td>Fastener specification database</td><td>Optional</td></tr>
       </table>
       
       <h3>Output Files</h3>
       <table>
       <tr><th>File Type</th><th>Extensions</th><th>Purpose</th><th>Content</th></tr>
       <tr><td>Modified BDF</td><td>.bdf</td><td>Updated model with new stiffness</td><td>Full FE model</td></tr>
       <tr><td>Updated Sets BDF</td><td>.bdf</td><td>Updated fastener group definitions</td><td>SET cards only</td></tr>
       <tr><td>Fastener Specs JSON</td><td>.json</td><td>Fastener specification export</td><td>Spec database</td></tr>
       </table>
       
       <h2>Nastran BDF Format</h2>
       
       <h3>Required Elements</h3>
       <h4>CBUSH Elements</h4>
       <div class="formula">
       CBUSH   EID     PID     GA      GB      X1      X2      X3      CID<br>
       +       G0      S       OCID    S1      S2      S3<br>
       <br>
       Example:<br>
       CBUSH   100001  1001    200001  200002  1.0     0.0     0.0
       </div>
       
       <h4>PBUSH Properties</h4>
       <div class="formula">
       PBUSH   PID     K1      K2      K3      K4      K5      K6      B1<br>
       +       B2      B3      B4      B5      B6      GE1     GE2     GE3<br>
       +       GE4     GE5     GE6     SA      ST      EA      ET      GS<br>
       <br>
       Example:<br>
       PBUSH   1001    1.0E7   1.0E7   1.0E7   100.0   1.0E8   1.0E8<br>
       </div>
       
       <h3>Supported Plate Elements</h3>
       <table>
       <tr><th>Element Type</th><th>Description</th><th>Nodes</th><th>Properties</th></tr>
       <tr><td>CQUAD4</td><td>4-node quadrilateral</td><td>4</td><td>PSHELL, PCOMP*</td></tr>
       <tr><td>CTRIA3</td><td>3-node triangle</td><td>3</td><td>PSHELL, PCOMP*</td></tr>
<tr><td>CQUAD8</td><td>8-node quadrilateral</td><td>8</td><td>PSHELL, PCOMP*</td></tr>
       <tr><td>CTRIA6</td><td>6-node triangle</td><td>6</td><td>PSHELL, PCOMP*</td></tr>
       </table>
       
       <h3>Property Cards</h3>
       <h4>PSHELL - Shell Properties</h4>
       <div class="formula">
       PSHELL  PID     MID1    T       MID2    12I/T³  MID3    TS/T    NSM<br>
       +       Z1      Z2      MID4<br>
       <br>
       Example:<br>
       PSHELL  2001    201     2.0     201     0.833   201     0.833
       </div>
       
       <h4>PCOMP - Composite Properties</h4>
       <div class="formula">
       PCOMP   PID     Z0      NSM     SB      FT      TREF    GE      LAM<br>
       +       MID1    T1      THETA1  SOUT1   MID2    T2      THETA2  SOUT2<br>
       +       ...continuing for each ply...<br>
       <br>
       Example:<br>
       PCOMP   3001    0.0                             0.0<br>
       +       301     0.125   0.0             301     0.125   45.0<br>
       +       301     0.125   -45.0           301     0.125   90.0<br>
       </div>
       
       <h4>PCOMPG - Global Ply Composite</h4>
       <div class="formula">
       PCOMPG  PID     Z0      NSM     SB      FT      TREF    GE      LAM<br>
       +       GPLYID1 GPLYID2 GPLYID3 ...<br>
       <br>
       References global ply definitions<br>
       </div>
       
       <h3>Material Cards</h3>
       <h4>MAT1 - Isotropic Material</h4>
       <div class="formula">
       MAT1    MID     E       G       NU      RHO     A       TREF    GE<br>
       +       ST      SC      SS      MCSID<br>
       <br>
       Example:<br>
       MAT1    201     70000.0 26900.0 0.3     2.7E-9<br>
       </div>
       
       <h4>MAT8 - Orthotropic Material</h4>
       <div class="formula">
       MAT8    MID     E1      E2      NU12    G12     G1Z     G2Z     RHO<br>
       +       A1      A2      TREF    XT      XC      YT      YC      S<br>
       +       GE      F12     STRN<br>
       <br>
       Example:<br>
       MAT8    301     150000. 9000.   0.3     4500.   4500.   3000.   1.6E-9<br>
       </div>
       
       <h2>SET Card Format</h2>
       
       <h3>Basic SET Definition</h3>
       <div class="formula">
       SET n = list of IDs<br>
       <br>
       Examples:<br>
       SET 100 = 1001, 1002, 1003, 1004<br>
       SET 101 = 2001 THRU 2010, 2015, 2020 THRU 2025<br>
       </div>
       
       <h3>SET with Naming (HyperMesh Style)</h3>
       <div class="formula">
       $HMSET setid type "setname" color<br>
       $HMSETTYPE setid "settype" color<br>
       <br>
       Example:<br>
       SET 100 = 1001, 1002, 1003<br>
       $HMSET      100        2 "WING_BOLTS__BOLT__6.35__NAS6204" 18<br>
       $HMSETTYPE  100 "regular" 18<br>
       </div>
       
       <h3>SET Naming Convention</h3>
       <div class="formula">
       Pattern: "DESCRIPTION1__DESCRIPTION2__DIAMETER__SPECIFICATION__ROWNUMBER__CPNUMBER__ENUM"<br>
       <br>
       Components:<br>
       - DESCRIPTION: General description of fastener group<br>
       - DIAMETER: Diameter in mm (e.g., 6.35)<br>
       - SPECIFICATION: Fastener spec (e.g., ELS438, MS20470)<br>
       - ROWNUMBER: Row number (e.g., 1, 2, 3, etc.)<br>
       - CPNUMBER: Connected plate number (e.g., 2, 3, 4, etc.)<br>
       - ENUM: Enumerator (e.g., 0, 1, 2 etc.)<br>
       <br>
       Examples:<br>
       "200__DuctRing.LHS__6.35__HST756__2__2__0"<br>
       "200__200Splice...WingUpperSkin.LHS__7.92__NAS1151V__1__3__0"
       </div>
       
       <h2>RBE3 Constraint Format</h2>
       
       <h3>RBE3 Definition</h3>
       <div class="formula">
       RBE3    EID     BLANK   REFGRID REFC    WT1     C1      G1      G2<br>
       +       G3      WT2     C2      G4      G5      G6      WT3     C3<br>
       +       G7      G8      G9      ...<br>
       <br>
       Purpose: Connect CBUSH nodes to plate elements
       </div>
       
       <h3>Typical CBUSH-RBE3-Plate Connection</h3>
       <div class="formula">
       Node Chain:
       CBUSH_Node → RBE3_Dependent → RBE3_Independent → Plate_Nodes<br>
       <br>
       Example:<br>
       CBUSH   100001  1001    200001  200002<br>
       RBE3    300001          200001  123     1.0     123     201001  201002<br>
       +       201003  201004<br>
       CQUAD4  400001  2001    201001  201002  201003  201004<br>
       </div>
       
       <h2>JSON Specification Format</h2>
       
       <h3>Fastener Specification Structure</h3>
       <div class="formula">
       {<br>
         "specification_name": {<br>
           "spec_name": "specification_name",<br>
           "fastener_type": "Bolt|Rivet|Hi-Lok|Blind Rivet",<br>
           "material": "Titanium|Steel|Aluminum|Custom",<br>
           "elastic_modulus": numeric_value_in_GPa<br>
         }<br>
       }
       </div>
       
       <h2>File Processing Workflow</h2>
       
       <h3>Input Processing</h3>
       <ol>
           <li><strong>BDF Parsing:</strong> PyNastran library reads main BDF file</li>
           <li><strong>Element Extraction:</strong> Identify CBUSH elements and properties</li>
           <li><strong>SET Parsing:</strong> Extract fastener groups from SET cards</li>
           <li><strong>Property Analysis:</strong> Extract plate and material properties</li>
           <li><strong>Connection Building:</strong> Identify multi-plate connections</li>
       </ol>
       
       <h3>Output Generation</h3>
       <ol>
           <li><strong>Property Update:</strong> Modify PBUSH stiffness values</li>
           <li><strong>Property Cloning:</strong> Create unique PIDs when needed (two CBUSH share the same PBUSH property)</li>
           <li><strong>BDF Writing:</strong> Generate updated BDF file</li>
           <li><strong>SET Export:</strong> Export updated fastener groups</li>
           <li><strong>Specification Save:</strong> Update specification database</li>
       </ol>
       
       <h2>Error Handling</h2>
       
       <h3>Common File Issues</h3>
       <table>
       <tr><th>Issue</th><th>Cause</th><th>Solution</th></tr>
       <tr><td>Parse errors</td><td>Invalid BDF syntax</td><td>Fix source file syntax</td></tr>
       <tr><td>Missing references</td><td>Undefined PIDs/MIDs</td><td>Complete property definitions</td></tr>
       <tr><td>Large file crashes</td><td>Insufficient memory</td><td>Increase RAM or reduce model</td></tr>
       <tr><td>Slow processing</td><td>Large model size</td><td>Use parallel processing</td></tr>
       </table>
       
       <div class="warning">
       <strong>Important:</strong> Always backup original files before processing. The tool creates modified copies, but backup ensures data safety in case of unexpected issues.
       </div>
       """
   
    elif topic_id == "coordinate_systems":
       return base_style + """
       <h1>Coordinate Systems and Transformations</h1>
       
       <h2>Overview</h2>
       <p>Understanding coordinate systems is crucial for accurate material orientation and fastener stiffness calculation. The tool handles multiple coordinate systems and their transformations automatically.</p>
       
       <h2>Coordinate System Hierarchy</h2>
       
       <h3>1. Global Coordinate System</h3>
       <ul>
           <li><strong>Definition:</strong> Model's primary coordinate system</li>
           <li><strong>Origin:</strong> (0, 0, 0) in model space</li>
           <li><strong>Axes:</strong> X, Y, Z as defined in BDF file</li>
           <li><strong>Usage:</strong> All node positions, element geometry</li>
       </ul>
       
       <h3>2. Element Local Coordinate System</h3>
       <ul>
           <li><strong>Definition:</strong> Local system for each plate element</li>
           <li><strong>Origin:</strong> Element centroid or first node</li>
           <li><strong>X-axis:</strong> Element edge direction (implementation dependent)</li>
           <li><strong>Z-axis:</strong> Normal to element surface</li>
       </ul>
       
       <h3>3. Material Coordinate System</h3>
       <ul>
           <li><strong>Definition:</strong> Principal material directions</li>
           <li><strong>Origin:</strong> Coincident with element system</li>
           <li><strong>Orientation:</strong> Defined by theta_mcid parameter</li>
           <li><strong>Usage:</strong> Material property directions (fiber directions)</li>
       </ul>
       
       <h3>4. Fastener Coordinate System</h3>
       <ul>
           <li><strong>Definition:</strong> Local system for CBUSH element</li>
           <li><strong>Origin:</strong> CBUSH centroid</li>
           <li><strong>e₁-axis:</strong> Along fastener (node1 to node2)</li>
           <li><strong>e₂, e₃-axes:</strong> Perpendicular shear directions</li>
       </ul>
       
       <h2>Element Coordinate System Determination</h2>
       
       <h3>CTRIA3 Elements</h3>
       <div class="formula">
       Nodes: G1, G2, G3 (corner nodes)<br>
       
       X-direction: normalize(G2 - G1)    // Side 1-2<br>
       Temp_Y: normalize(G3 - G1)         // Side 1-3  <br>
       Z-direction: normalize(X × Temp_Y) // Normal (right-hand rule)<br>
       Y-direction: Z × X                 // Complete right-hand system<br>
       </div>
       
       <h3>CQUAD4 Elements</h3>
       <div class="formula">
       Nodes: G1, G2, G3, G4 (corner nodes in order)<br>
       
       X-direction: normalize(G2 - G1)    // Side 1-2 (simplified)<br>
       Temp_Y: normalize(G4 - G1)         // Side 1-4<br>
       Z-direction: normalize(X × Temp_Y) // Normal (right-hand rule)<br>
       Y-direction: Z × X                 // Complete right-hand system<br>
       <br>
       Note: Full implementation should use diagonal bisection method
       </div>
       
       <h3>Higher-Order Elements</h3>
       <ul>
           <li><strong>CTRIA6:</strong> Use first 3 corner nodes (ignore mid-side nodes)</li>
           <li><strong>CQUAD8:</strong> Use first 4 corner nodes (ignore mid-side nodes)</li>
           <li><strong>Approach:</strong> Degenerate to linear element for coordinate calculation</li>
       </ul>
       
       <h2>Material Orientation (theta_mcid)</h2>
       
       <h3>Direct Angle Specification</h3>
       <div class="formula">
       theta_mcid = angle_in_degrees (floating point)<br>
       <br>
       Rotation: Rotate element X-axis by theta about element Z-axis<br>
       Material_X = rotate(Element_X, Element_Z, theta_radians)<br>
       <br>
       Positive angle: Counter-clockwise rotation when looking down element Z-axis<br>
       </div>
       
       <h3>Coordinate System Reference</h3>
       <div class="formula">
       theta_mcid = coordinate_system_ID (positive integer)<br>
       <br>
       Process:<br>
       1. Look up coordinate system definition<br>
       2. Extract X-axis direction in global coordinates<br>
       3. Project onto element surface (remove Z-component)<br>
       4. Normalize to get material direction<br>
       </div>
       
       <h3>Zero or Undefined Orientation</h3>
       <div class="formula">
       theta_mcid = 0.0 or undefined<br>
       <br>
       Result: Material X-direction = Element X-direction (no rotation)
       </div>
       
       <h2>Fastener Coordinate System</h2>
       
       <h3>Primary Axis (e₁)</h3>
       <div class="formula">
       e₁ = normalize(Node2_position - Node1_position)<br>
       <br>
       This is the fastener axial direction (for K₁ stiffness)
       </div>
       
       <h3>Orientation Vector Determination</h3>
       <p>The tool determines fastener orientation using available CBUSH orientation data:</p>
       
       <h4>Method 1: X Vector (Most Common)</h4>
       <div class="formula">
       IF CBUSH.x is defined:<br>
       >    orientation_vector = CBUSH.x<br>
       </div>
       
       <h4>Method 2: G0 Node</h4>
       <div class="formula">
       IF CBUSH.g0 is defined:<br>
       >    orientation_vector = Node_G0_position - Node1_position
       </div>
       
       <h4>Method 3: Default Orientation</h4>
       <div class="formula">
       IF no orientation specified:<br>
       >    Use default perpendicular vector to e₁<br>
       >    Select vector with minimum component overlap
       </div>
       
       <h3>Secondary Axes (e₂, e₃)</h3>
       <div class="formula">
       e₃ = normalize(e₁ × orientation_vector)<br>
       e₂ = e₃ × e₁<br>
       <br>
       This creates a right-handed coordinate system:<br>
       - e₁: Axial direction (tension/compression)<br>
       - e₂: Shear direction 1 (K₂ stiffness)<br>
       - e₃: Shear direction 2 (K₃ stiffness)<br>
       </div>
       
       <h2>Coordinate Transformation Process</h2>
       
       <h3>Step 1: Element Geometry Analysis</h3>
       <div class="formula">
       # Get element nodes and calculate local system<br>
       element_nodes = get_element_nodes(element_id)<br>
       element_x, element_y, element_z = calculate_element_system(element_nodes)<br>
       </div>
       
       <h3>Step 2: Material Direction Calculation</h3>
       <div class="formula">
       # Apply material orientation<br>
       theta = get_theta_mcid(element)<br>
       material_x = apply_material_rotation(element_x, element_z, theta)<br>
       </div>
       
       <h3>Step 3: Fastener System Construction</h3>
       <div class="formula">
       # Build fastener coordinate system<br>
       fastener_e1, fastener_e2, fastener_e3 = get_fastener_system(cbush_element)<br>
       </div>
       
       <h3>Step 4: Material-Fastener Angle Calculation</h3>
       <div class="formula">
       # Project fastener directions onto element surface<br>
       fastener_e2_projected = project_onto_plane(fastener_e2, element_z)<br>
       fastener_e3_projected = project_onto_plane(fastener_e3, element_z)<br>
       <br>
       # Calculate angles<br>
       angle_e2 = angle_between_vectors(material_x, fastener_e2_projected)<br>
       angle_e3 = angle_between_vectors(material_x, fastener_e3_projected)<br>
       </div>
       
       <h2>Transformation Utilities</h2>
       
       <h3>Vector Projection onto Plane</h3>
       <div class="formula">
       def project_onto_plane(vector, normal):<br>
           # Remove component parallel to normal<br>
           normal = normalize(normal)<br>
           parallel_component = dot(vector, normal) * normal<br>
           projected = vector - parallel_component<br>
           return normalize(projected)
       </div>
       
       <h3>Rotation in Plane</h3>
       <div class="formula">
       def rotate_in_plane(vector, normal, angle):<br>
           # Rodrigues' rotation formula<br>
           cos_a = cos(angle)<br>
           sin_a = sin(angle)<br>
           cross_prod = cross(normal, vector)<br>
           dot_prod = dot(normal, vector)<br>
           <br>
           rotated = (vector * cos_a + <br>
                     cross_prod * sin_a + <br>
                     normal * dot_prod * (1 - cos_a))<br>
           return normalize(rotated)
       </div>
       
       <h3>Angle Between Vectors</h3>
       <div class="formula">
       def angle_between_vectors(v1, v2):<br>
           v1_norm = normalize(v1)<br>
           v2_norm = normalize(v2)<br>
           cos_angle = clip(dot(v1_norm, v2_norm), -1.0, 1.0)<br>
           return arccos(abs(cos_angle))  # Always positive angle
       </div>

       """
   
#     elif topic_id == "clt_theory":
#        return base_style + """
#        <h1>Classical Lamination Theory Implementation</h1>
       
#        <h2>Overview</h2>
#        <p>Classical Lamination Theory (CLT) is implemented to calculate accurate effective properties for composite laminates. This is essential for proper fastener stiffness calculation in composite structures.</p>
       
#        <h2>Laminate Property Cards</h2>
       
#        <h3>PCOMP - Ply-by-Ply Definition</h3>
#        <div class="formula">
#        Structure:
#        - Layer-by-layer material and angle specification
#        - Individual ply thicknesses
#        - Symmetry options (SYM, SYMEM, etc.)
#        - Direct material reference per ply
#        </div>
       
#        <h3>PCOMPG - Global Ply Reference</h3>
#        <div class="formula">
#        Structure:
#        - References to global ply database
#        - Ply IDs instead of material IDs
#        - Shared ply definitions across elements
#        - More complex data structure
#        </div>
       
#        <h3>PCOMPP - Parametric Definition</h3>
#        <div class="formula">
#        Structure:
#        - Parametric composite definition
#        - References to base PCOMP properties
#        - Scaling factors and modifications
#        - Less commonly used
#        </div>
       
#        <h2>CLT Mathematical Foundation</h2>
       
#        <h3>Ply Stiffness Matrix (Q)</h3>
#        <div class="formula">
#        For orthotropic ply (material coordinates):
       
#        Q₁₁ = E₁ / (1 - ν₁₂ν₂₁)
#        Q₂₂ = E₂ / (1 - ν₁₂ν₂₁)  
#        Q₁₂ = ν₁₂E₂ / (1 - ν₁₂ν₂₁)
#        Q₆₆ = G₁₂
#        Q₁₆ = Q₂₆ = 0 (for orthotropic materials)
       
#        Where: ν₂₁ = ν₁₂E₂/E₁
#        </div>
       
#        <h3>Transformed Stiffness Matrix (Q̄)</h3>
#        <div class="formula">
#        For ply rotated by angle θ:
       
#        Q̄₁₁ = Q₁₁c⁴ + 2(Q₁₂ + 2Q₆₆)c²s² + Q₂₂s⁴
#        Q̄₂₂ = Q₁₁s⁴ + 2(Q₁₂ + 2Q₆₆)c²s² + Q₂₂c⁴
#        Q̄₁₂ = (Q₁₁ + Q₂₂ - 4Q₆₆)c²s² + Q₁₂(c⁴ + s⁴)
#        Q̄₁₆ = (Q₁₁ - Q₁₂ - 2Q₆₆)c³s + (Q₁₂ - Q₂₂ + 2Q₆₆)cs³
#        Q̄₂₆ = (Q₁₁ - Q₁₂ - 2Q₆₆)cs³ + (Q₁₂ - Q₂₂ + 2Q₆₆)c³s
#        Q̄₆₆ = (Q₁₁ + Q₂₂ - 2Q₁₂ - 2Q₆₆)c²s² + Q₆₆(c⁴ + s⁴)
       
#        Where: c = cos(θ), s = sin(θ)
#        </div>
       
#        <h3>ABD Matrix Construction</h3>
#        <div class="formula">
#        Aᵢⱼ = Σ(Q̄ᵢⱼ)ₖ(zₖ - zₖ₋₁)           [Extensional stiffness]
#        Bᵢⱼ = ½Σ(Q̄ᵢⱼ)ₖ(zₖ² - zₖ₋₁²)         [Coupling stiffness]  
#        Dᵢⱼ = ⅓Σ(Q̄ᵢⱼ)ₖ(zₖ³ - zₖ₋₁³)         [Bending stiffness]
       
#        Where:
#        k = ply number
#        z = distance from laminate midplane
#        zₖ₋₁ = bottom of ply k
#        zₖ = top of ply k
#        </div>
       
#        <h2>Implementation Details</h2>
       
#        <h3>Ply Extraction Process</h3>
#        <h4>PCOMP Processing</h4>
#        <div class="formula">
#        for i in range(n_plies):
#            ply = {
#                'material_id': prop.mids[i],
#                'thickness': prop.thicknesses[i], 
#                'angle': prop.thetas[i],
#                'ply_id': i + 1
#            }
#            plies.append(ply)
       
#        # Handle symmetry
#        if prop.lam in ['SYM', 'SYMEM']:
#            symmetric_plies = create_symmetric_plies(plies)
#            plies.extend(symmetric_plies)
#        </div>
       
#        <h4>PCOMPG Processing</h4>
#        <div class="formula">
#        # Extract from global ply database
#        for ply_id in element.plylist:
#            global_ply = model.plies[ply_id]
#            ply = {
#                'material_id': global_ply.mid,
#                'thickness': global_ply.t,
#                'angle': global_ply.theta,
#                'ply_id': ply_id,
#                'global_ply': True
#            }
#            plies.append(ply)
#        </div>
       
#        <h3>Z-Coordinate Calculation</h3>
#        <div class="formula">
#        total_thickness = sum(ply['thickness'] for ply in plies)
#        z_bottom = -total_thickness / 2.0  # Start from bottom
       
#        for ply in plies:
#            ply['z_bottom'] = z_bottom
#            ply['z_top'] = z_bottom + ply['thickness']
#            z_bottom = ply['z_top']
#        </div>
       
#        <h3>Material Property Lookup</h3>
#        <div class="formula">
#        def get_material_properties(material_id):
#            material = model.materials[material_id]
           
#            if material.type == 'MAT8':
#                return {
#                    'E1': material.e11,
#                    'E2': material.e22, 
#                    'G12': material.g12,
#                    'NU12': material.nu12
#                }
#            elif material.type == 'MAT1':
#                return {
#                    'E1': material.e,
#                    'E2': material.e,
#                    'G12': material.g, 
#                    'NU12': material.nu
#                }
#        </div>
       
#        <h2>Effective Properties Calculation</h2>
       
#        <h3>In-Plane Properties from A Matrix</h3>
#        <div class="formula">
#        # Invert A matrix to get compliance
#        a = inverse(A)
       
#        # Calculate effective moduli
#        Ex = 1 / (a₁₁ × total_thickness)
# Ey = 1 / (a₂₂ × total_thickness)
#        Gxy = 1 / (a₆₆ × total_thickness)
       
#        # Effective Poisson's ratio
#        νxy = -a₁₂ / a₁₁
       
#        # Geometric mean for general use
#        E_effective = √(Ex × Ey)
#        </div>
       
#        <h3>Weighted Average Fallback</h3>
#        <div class="formula">
#        # If A matrix inversion fails, use weighted average
#        Ex_avg = Σ(Ex_ply × t_ply) / total_thickness
#        Ey_avg = Σ(Ey_ply × t_ply) / total_thickness  
#        Gxy_avg = Σ(Gxy_ply × t_ply) / total_thickness
       
#        Where ply properties use Jones transformation:
#        Ex_ply = jones_transform(E1, E2, G12, ν12, θ)
#        </div>
       
#        <h2>Symmetry Handling</h2>
       
#        <h3>Symmetric Laminates</h3>
#        <div class="formula">
#        LAM = 'SYM':  Full symmetry about midplane
#        LAM = 'SYMEM': Symmetric but exclude midplane ply
#        LAM = 'SYBEND': Symmetric for bending only
#        LAM = 'SYSMEAR': Symmetric with smeared properties
       
#        Implementation:
#        if is_symmetric:
#            for ply in reversed(base_plies):
#                symmetric_ply = ply.copy()
#                symmetric_ply['ply_id'] = next_ply_id
#                symmetric_plies.append(symmetric_ply)
#        </div>
       
#        <h3>Asymmetric Laminates</h3>
#        <ul>
#            <li>No automatic ply duplication</li>
#            <li>Use exactly as specified in PCOMP</li>
#            <li>May have coupling (B matrix non-zero)</li>
#            <li>Requires careful handling in analysis</li>
#        </ul>
       
#        <h2>Coupling Effects</h2>
       
#        <h3>Extension-Bending Coupling</h3>
#        <div class="formula">
#        If B ≠ 0: Extension and bending are coupled
       
#        Physical meaning:
#        - In-plane loads cause out-of-plane deformation
#        - Bending moments cause in-plane strains
#        - Common in asymmetric laminates
       
#        For fastener analysis:
#        - Usually neglect coupling effects
#        - Focus on in-plane properties (A matrix)
#        - Use effective in-plane moduli
#        </div>
       
#        <h3>Coupling Indicators</h3>
#        <div class="formula">
#        Coupling ratio = max(|Bᵢⱼ|) / max(|Aᵢⱼ| × thickness)
       
#        Typical guidelines:
#        - Ratio < 0.1: Weakly coupled, safe to ignore
#        - Ratio 0.1-0.3: Moderately coupled, use caution
#        - Ratio > 0.3: Strongly coupled, detailed analysis needed
#        </div>
       
#        <h2>Special Laminate Cases</h2>
       
#        <h3>Quasi-Isotropic Laminates</h3>
#        <div class="formula">
#        Common stacking: [0/±45/90]ₛ or similar
       
#        Properties:
#        - Ex ≈ Ey (nearly isotropic in-plane)
#        - Minimal coupling if symmetric
#        - Good general-purpose properties
       
#        CLT validation:
#        - Check Ex/Ey ratio ≈ 1.0
#        - Verify A₁₆ ≈ A₂₆ ≈ 0
#        </div>
       
#        <h3>Unidirectional Laminates</h3>
#        <div class="formula">
#        Single ply angle: [0]ₙ or [90]ₙ
       
#        Properties:
#        - Highly orthotropic (Ex >> Ey or vice versa)
#        - Orientation critical for stiffness
#        - Sensitive to fiber direction accuracy
       
#        Special handling:
#        - Use actual Ex, Ey (not geometric mean)
#        - Material orientation becomes critical
#        </div>
       
#        <h3>Fabric Laminates</h3>
#        <div class="formula">
#        Woven or braided fabrics: Balanced properties
       
#        Modeling approaches:
#        - Single equivalent layer with smeared properties
#        - Individual ply modeling of warp/weft
#        - Use MAT8 with balanced E1 ≈ E2
#        </div>
       
#        <h2>Validation and Quality Checks</h2>
       
#        <h3>Physical Reasonableness</h3>
#        <div class="formula">
#        Checks to perform:
       
#        1. Positive definite A matrix
#        2. Ex, Ey > 0
#        3. |νxy| < 1.0
#        4. Gxy > 0
#        5. Reasonable property ratios:
#           - 0.1 < Ex/Ey < 10 (for most laminates)
#           - Gxy < min(Ex, Ey)
#        </div>
       
#        <h3>CLT Matrix Validation</h3>
#        <div class="formula">
#        Matrix properties:
#        - A matrix: 3×3, symmetric, positive definite
#        - B matrix: 3×3, symmetric (may be zero)
#        - D matrix: 3×3, symmetric, positive definite
       
#        Symmetry check: |Aᵢⱼ - Aⱼᵢ| < tolerance
#        </div>
       
#        <h3>Ply-Level Validation</h3>
#        <ul>
#            <li><strong>Thickness:</strong> All ply thicknesses > 0</li>
#            <li><strong>Materials:</strong> Valid material references</li>
#            <li><strong>Angles:</strong> Reasonable ply angles (-90° to +90°)</li>
#            <li><strong>Stacking:</strong> Proper z-coordinate progression</li>
#        </ul>
       
#        <h2>Error Handling and Fallbacks</h2>
       
#        <h3>Matrix Inversion Failures</h3>
#        <div class="formula">
#        Primary method: A_compliance = inverse(A_stiffness)
       
#        If inversion fails:
#        1. Check matrix conditioning
#        2. Try pseudo-inverse
#        3. Fall back to weighted average method
#        4. Use simplified isotropic equivalent
#        </div>
       
#        <h3>Missing Data Handling</h3>
#        <ul>
#            <li><strong>Missing materials:</strong> Use default MAT8 properties</li>
#            <li><strong>Missing plies:</strong> Skip and warn user</li>
#            <li><strong>Invalid angles:</strong> Use 0° default</li>
#            <li><strong>Zero thickness:</strong> Skip ply in calculations</li>
#        </ul>
       
#        <h3>Numerical Issues</h3>
#        <div class="formula">
#        Common problems and solutions:
       
#        1. Very thin plies: tₘᵢₙ = 1e-6 mm
#        2. Large angle variations: Use consistent units (degrees)
#        3. Extreme property ratios: E₁/E₂ > 1000 → numerical issues
#        4. Nearly singular matrices: Check condition number
#        </div>
       
#        <h2>Integration with Fastener Analysis</h2>
       
#        <h3>Effective Properties for Huth Formula</h3>
#        <div class="formula">
#        From CLT analysis:
#        Ex, Ey, Gxy, νxy, total_thickness
       
#        For material orientation approach:
#        - Use Ex, Ey directly with Jones transformation
#        - Apply orientation relative to fastener coordinates
       
#        For effective modulus approach:
#        - Use E_effective = √(Ex × Ey)
#        - Single value for both K₂ and K₃
#        </div>
       
#        <h3>Composite Material Parameters</h3>
#        <div class="formula">
#        Huth formula modifications for composites:
       
#        b_composite = 4.2 (vs 3.0 for metals)
       
#        Accounts for:
#        - Different bearing behavior
#        - Anisotropic material response
#        - Fiber-matrix interaction effects
#        </div>
       
#        <h2>Advanced Features</h2>
       
#        <h3>Temperature Effects</h3>
#        <ul>
#            <li>Material property temperature dependence</li>
#            <li>Thermal expansion mismatch effects</li>
#            <li>Temperature-dependent material degradation</li>
#            <li>Currently not implemented (future enhancement)</li>
#        </ul>
       
#        <h3>Moisture Effects</h3>
#        <ul>
#            <li>Hygroscopic swelling in composites</li>
#            <li>Property degradation with moisture</li>
#            <li>Environmental conditioning effects</li>
#            <li>Currently not implemented (future enhancement)</li>
#        </ul>
       
#        <h3>Progressive Damage</h3>
#        <ul>
#            <li>Ply-by-ply failure analysis</li>
#            <li>Property reduction with damage</li>
#            <li>Matrix cracking and delamination</li>
#            <li>Beyond current scope (specialized analysis)</li>
#        </ul>
       
#        <div class="note">
#        <strong>Implementation Note:</strong> The CLT implementation focuses on linear elastic analysis for fastener stiffness calculation. Advanced effects like damage, temperature, and moisture are not currently included but could be added as enhancements.
#        </div>
#        """
   
    # elif topic_id == "validation":
    #    return base_style + """
    #    <h1>Validation and Quality Assurance</h1>
       
    #    <h2>Overview</h2>
    #    <p>Comprehensive validation ensures the accuracy and reliability of fastener stiffness calculations. This section outlines validation methods, quality checks, and verification procedures.</p>
       
    #    <h2>Multi-Level Validation Strategy</h2>
       
    #    <h3>Level 1: Input Validation</h3>
    #    <ul>
    #        <li><strong>File Format:</strong> Valid BDF syntax and structure</li>
    #        <li><strong>Element Completeness:</strong> All required cards present</li>
    #        <li><strong>Reference Integrity:</strong> Valid node, property, and material references</li>
    #        <li><strong>Geometric Validity:</strong> Non-degenerate elements and reasonable dimensions</li>
    #    </ul>
       
    #    <h3>Level 2: Connection Validation</h3>
    #    <ul>
    #        <li><strong>CBUSH-Plate Connectivity:</strong> Proper RBE3 connections</li>
    #        <li><strong>Multi-Plate Identification:</strong> Correct grouping of connected fasteners</li>
    #        <li><strong>Plate Ordering:</strong> Logical sequence through thickness</li>
    #        <li><strong>Property Extraction:</strong> Successful material property retrieval</li>
    #    </ul>
       
    #    <h3>Level 3: Calculation Validation</h3>
    #    <ul>
    #        <li><strong>Stiffness Magnitude:</strong> Physically reasonable values</li>
    #        <li><strong>Consistency Checks:</strong> Similar fasteners have similar stiffness</li>
    #        <li><strong>Formula Verification:</strong> Huth formula implementation correctness</li>
    #        <li><strong>Coordinate Transformation:</strong> Proper orientation handling</li>
    #    </ul>
       
    #    <h3>Level 4: Output Validation</h3>
    #    <ul>
    #        <li><strong>Model Integrity:</strong> Updated model remains valid</li>
    #        <li><strong>Property Uniqueness:</strong> No unintended property sharing</li>
    #        <li><strong>Completeness:</strong> All target elements updated</li>
    #        <li><strong>Documentation:</strong> Proper logging and traceability</li>
    #    </ul>
       
    #    <h2>Automated Quality Checks</h2>
       
    #    <h3>Stiffness Value Validation</h3>
    #    <div class="formula">
    #    Typical stiffness ranges for aerospace fasteners:
       
    #    K₁ (Axial):     1E6 - 1E8 N/mm
    #    K₂, K₃ (Shear): 1E5 - 1E7 N/mm  
    #    K₄ (Torsion):   10 - 1000 N⋅mm/rad
    #    K₅, K₆ (Bend):  1E7 - 1E9 N⋅mm/rad
       
    #    Warnings issued for values outside these ranges
    #    </div>
       
    #    <h3>Consistency Checks</h3>
    #    <div class="formula">
    #    For fasteners with identical:
    #    - Diameter
    #    - Material specification  
    #    - Plate materials and thicknesses
    #    - Connection type
       
    #    Stiffness values should be within ±5% tolerance
    #    </div>
       
    #    <h3>Physical Relationship Validation</h3>
    #    <div class="formula">
    #    Expected relationships:
       
    #    1. K₁ > K₂, K₃ (typically by factor of 2-10)
    #    2. K₅, K₆ >> K₂, K₃ (bending much stiffer than shear)
    #    3. Larger diameter → higher stiffness
    #    4. Higher modulus materials → higher stiffness
    #    5. Thicker plates → lower fastener stiffness
    #    </div>
       
    #    <h2>Test Case Validation</h2>
       
    #    <h3>Simple Two-Plate Test Cases</h3>
    #    <h4>Test Case 1: Identical Aluminum Plates</h4>
    #    <div class="formula">
    #    Configuration:
    #    - Plate 1: 2.0mm aluminum, E=70000 MPa
    #    - Plate 2: 2.0mm aluminum, E=70000 MPa  
    #    - Fastener: 6.35mm titanium bolt, E=110000 MPa
       
    #    Expected results:
    #    - K₁ ≈ 3.47E7 N/mm
    #    - K₂ = K₃ ≈ 4.2E6 N/mm (symmetric case)
    #    - Material orientation should not affect isotropic case
    #    </div>
       
    #    <h4>Test Case 2: Dissimilar Materials</h4>
    #    <div class="formula">
    #    Configuration:
    #    - Plate 1: 1.5mm steel, E=210000 MPa
    #    - Plate 2: 3.0mm aluminum, E=70000 MPa
    #    - Fastener: 6.35mm steel bolt, E=210000 MPa
       
    #    Expected results:
    #    - Higher overall stiffness due to steel
    #    - Asymmetric response due to material mismatch
    #    - Validation against hand calculations
    #    </div>
       
    #    <h3>Composite Material Test Cases</h3>
    #    <h4>Test Case 3: Unidirectional Composite</h4>
    #    <div class="formula">
    #    Configuration:
    #    - Single ply: [0°] IM7/8552, t=0.125mm
    #    - Material: E₁=161000, E₂=11380, G₁₂=5170, ν₁₂=0.32
    #    - Fastener aligned with fibers (0°)
       
    #    Expected results:
    #    - High stiffness in fiber direction
    #    - Lower stiffness perpendicular to fibers
    #    - Strong orientation dependence
    #    </div>
       
    #    <h4>Test Case 4: Quasi-Isotropic Laminate</h4>
    #    <div class="formula">
    #    Configuration:
    #    - Stacking: [0/45/-45/90]ₛ IM7/8552
    #    - Total thickness: 1.0mm
    #    - Fastener at various orientations
       
    #    Expected results:
    #    - Nearly isotropic in-plane properties
    #    - Minimal orientation dependence
    #    - Ex ≈ Ey ≈ 55000 MPa (typical)
    #    </div>
       
    #    <h2>Comparison with Reference Solutions</h2>
       
    #    <h3>Analytical Solutions</h3>
    #    <ul>
    #        <li><strong>Simple configurations:</strong> Hand calculations using Huth formula</li>
    #        <li><strong>Isotropic cases:</strong> Closed-form solutions available</li>
    #        <li><strong>Symmetric laminates:</strong> CLT-based calculations</li>
    #        <li><strong>Published data:</strong> Aerospace handbook values</li>
    #    </ul>
       
    #    <h3>Commercial Software Comparison</h3>
    #    <table>
    #    <tr><th>Software</th><th>Capabilities</th><th>Comparison Method</th></tr>
    #    <tr><td>HyperMesh</td><td>Fastener modeling</td><td>Direct stiffness comparison</td></tr>
    #    <tr><td>Patran</td><td>Multi-point constraints</td><td>Equivalent spring properties</td></tr>
    #    <tr><td>FEMAP</td><td>Connection elements</td><td>Modal frequency comparison</td></tr>
    #    <tr><td>Custom tools</td><td>Company-specific</td><td>Legacy model validation</td></tr>
    #    </table>
       
    #    <h3>Experimental Validation</h3>
    #    <ul>
    #        <li><strong>Fastener pull-out tests:</strong> Axial stiffness validation</li>
    #        <li><strong>Lap shear tests:</strong> Shear stiffness validation</li>
    #        <li><strong>Modal testing:</strong> Overall structural response</li>
    #        <li><strong>Correlation studies:</strong> Test-analysis correlation</li>
    #    </ul>
       
    #    <h2>Error Detection and Reporting</h2>
       
    #    <h3>Common Error Categories</h3>
    #    <table>
    #    <tr><th>Error Type</th><th>Detection Method</th><th>Resolution</th></tr>
    #    <tr><td>Missing properties</td><td>Reference checking</td><td>Use defaults with warning</td></tr>
    #    <tr><td>Degenerate elements</td><td>Geometric validation</td><td>Skip with error message</td></tr>
    #    <tr><td>Invalid materials</td><td>Property range checking</td><td>Use fallback properties</td></tr>
    #    <tr><td>Connection failures</td><td>Connectivity analysis</td><td>Manual review required</td></tr>
    #    </table>
       
    #    <h3>Warning System</h3>
    #    <div class="formula">
    #    Warning levels:
       
    #    INFO:    Informational messages, normal operation
    #    WARNING: Potential issues, results may be affected
    #    ERROR:   Serious problems, manual intervention needed
    #    CRITICAL: Fatal errors, processing cannot continue
    #    </div>
       
    #    <h3>Quality Metrics</h3>
    #    <ul>
    #        <li><strong>Success Rate:</strong> Percentage of elements successfully processed</li>
    #        <li><strong>Consistency Index:</strong> Variation in similar fastener stiffness</li>
    #        <li><strong>Coverage Ratio:</strong> Percentage of fasteners with specifications</li>
    #        <li><strong>Validation Score:</strong> Overall model quality assessment</li>
    #    </ul>
       
    #    <h2>Verification Procedures</h2>
       
    #    <h3>Pre-Processing Verification</h3>
    #    <ol>
    #        <li><strong>Model Review:</strong> Visual inspection of fastener locations</li>
    #        <li><strong>Connectivity Check:</strong> Verify CBUSH-plate connections</li>
    #        <li><strong>Material Audit:</strong> Review material property assignments</li>
    #        <li><strong>Specification Validation:</strong> Confirm fastener specifications</li>
    #    </ol>
       
    #    <h3>Post-Processing Verification</h3>
    #    <ol>
    #        <li><strong>Stiffness Review:</strong> Check calculated values for reasonableness</li>
    #        <li><strong>Consistency Analysis:</strong> Compare similar fasteners</li>
    #        <li><strong>Spot Checks:</strong> Manual verification of sample calculations</li>
    #        <li><strong>Model Validation:</strong> Ensure updated model loads correctly</li>
    #    </ol>
       
    #    <h3>Documentation Requirements</h3>
    #    <ul>
    #        <li><strong>Calculation Log:</strong> Complete record of all operations</li>
    #        <li><strong>Warning Summary:</strong> List of all warnings and errors</li>
    #        <li><strong>Quality Report:</strong> Statistical summary of results</li>
    #        <li><strong>Traceability Matrix:</strong> Element-to-calculation mapping</li>
    #    </ul>
       
    #    <h2>Best Practices for Validation</h2>
       
    #    <h3>Model Preparation</h3>
    #    <ul>
    #        <li>Start with clean, well-structured BDF files</li>
    #        <li>Use consistent naming conventions</li>
    #        <li>Validate input model before processing</li>
    #        <li>Backup original files before modification</li>
    #    </ul>
       
    #    <h3>Incremental Validation</h3>
    #    <ul>
    #        <li>Test on small subsets first</li>
    #        <li>Validate each processing stage</li>
    #        <li>Compare with previous results when available</li>
    #        <li>Document all assumptions and decisions</li>
    #    </ul>
       
    #    <h3>Quality Assurance Process</h3>
    #    <ol>
    #        <li><strong>Independent Review:</strong> Second engineer reviews setup</li>
    #        <li><strong>Sample Verification:</strong> Hand-check representative cases</li>
    #        <li><strong>Trend Analysis:</strong> Look for patterns in results</li>
    #        <li><strong>Sensitivity Study:</strong> Assess impact of parameter variations</li>
    #    </ol>
       
    #    <h2>Validation Metrics and Criteria</h2>
       
    #    <h3>Acceptance Criteria</h3>
    #    <div class="formula">
    #    Model acceptance requires:
       
    #    1. Success rate > 95%
    #    2. Stiffness values within expected ranges
    #    3. Consistency index < 10% for similar fasteners
    #    4. No critical errors
    #    5. All warnings reviewed and justified
    #    </div>
       
    #    <h3>Performance Metrics</h3>
    #    <table>
    #    <tr><th>Metric</th><th>Good</th><th>Acceptable</th><th>Needs Review</th></tr>
    #    <tr><td>Success Rate</td><td>>98%</td><td>95-98%</td><td><95%</td></tr>
    #    <tr><td>Consistency</td><td><5%</td><td>5-10%</td><td>>10%</td></tr>
    #    <tr><td>Coverage</td><td>>90%</td><td>80-90%</td><td><80%</td></tr>
    #    <tr><td>Warnings</td><td><10</td><td>10-50</td><td>>50</td></tr>
    #    </table>
       
    #    <div class="warning">
    #    <strong>Critical:</strong> Never use calculated stiffness values in critical analyses without proper validation. When in doubt, perform additional verification using independent methods or simplified conservative estimates.
    #    </div>
    #    """
   
#     elif topic_id == "troubleshooting":
#        return base_style + """
#        <h1>Troubleshooting Guide</h1>
       
#        <h2>Common Issues and Solutions</h2>
       
#        <h3>File Loading Problems</h3>
#        <table>
#        <tr><th>Error</th><th>Cause</th><th>Solution</th></tr>
#        <tr><td>"Invalid BDF format"</td><td>Syntax errors in BDF file</td><td>Check line formatting, field alignment</td></tr>
#        <tr><td>"Unknown card type"</td><td>Unsupported Nastran cards</td><td>Remove or comment unsupported cards</td></tr>
#        <tr><td>"Missing continuation"</td><td>Incomplete multi-line cards</td><td>Verify + continuation markers</td></tr>
#        <tr><td>"Memory error"</td><td>File too large for available RAM</td><td>Increase RAM or split model</td></tr>
#        </table>
       
#        <h3>Connection Identification Issues</h3>
       
#        <h4>No Connections Found</h4>
#        <div class="formula">
#        Possible causes:<br>
#        <br>
#        1. CBUSH elements not connected through shared nodes<br>
#        2. Missing RBE3 elements linking CBUSH to plates<br>
#        <br>
#        Debugging steps:<br>
#        1. Check CBUSH node IDs exist in model<br>
#        2. Verify RBE3 connections
#        </div>
       
#        <h4>Incorrect Multi-Plate Grouping</h4>
#        <div class="formula">
#        Symptoms:<br>
#        - Too many connections (should be fewer for multi-plate)<br>
#        - Single CBUSHes not grouped with neighbors<br>
#        - Inconsistent plate counts<br>
#        <br>
#        Solutions:<br>
#        1. Verify shared node connectivity<br>
#        2. Check RBE3 independent node assignments<br>
#        3. Review plate element connectivity<br>
#        4. Manual verification of sample connections<br>
#        </div>
       
#        <h3>Material Property Issues</h3>
       
#        <h4>Missing Material Properties</h4>
#        <table>
#        <tr><th>Issue</th><th>Detection</th><th>Solution</th></tr>
#        <tr><td>Undefined MID</td><td>Reference checking</td><td>Define missing materials or use defaults</td></tr>
#        <tr><td>Zero modulus</td><td>Property validation</td><td>Correct material definition</td></tr>
#        <tr><td>Invalid MAT8</td><td>Property extraction</td><td>Check E1, E2, G12, NU12 values</td></tr>
#        <tr><td>PCOMP errors</td><td>Ply extraction</td><td>Verify ply definitions and materials</td></tr>
#        </table>
       
#        <h4>Composite Property Problems</h4>
#        <div class="formula">
#        Common PCOMP issues:<br>
#        <br>
#        1. Missing ply materials → Use backup MAT8 properties<br>
#        2. Invalid ply angles → Default to 0°<br>
#        3. Zero ply thickness → Skip ply with warning<br>
#        4. Symmetry errors → Manual ply definition<br>
#        5. CLT matrix singular → Use weighted average fallback<br>
#        </div>
       
#        <h3>Calculation Errors</h3>
       
#        <h4>Unreasonable Stiffness Values</h4>
#        <div class="formula">
#        Value ranges for diagnosis:<br>
#        <br>
#        Too high (>1E9):<br>
#        - Check thickness units (should be mm)<br>
#        - Verify diameter units (should be mm)<br>
#        - Review material modulus units (should be MPa)<br>
#        <br>
#        Too low (<1E3):<br>
#        - Check for very thick plates<br>
#        - Verify fastener diameter not too small<br>
#        - Review material properties for errors<br>
#        <br>
#        Negative values:<br>
#        - Indicates calculation error<br>
#        - Check input parameters<br>
#        - Review coordinate transformations
#        </div>
       
#        <h4>Material Orientation Problems</h4>
#        <div class="formula">
#        Symptoms:<br>
#        - K2 ≠ K3 when they should be equal<br>
#        - Unexpected orientation sensitivity<br>
#        - Inconsistent results for similar elements<br>
#        <br>
#        Debugging:<br>
#        1. Check theta_mcid values on elements<br>
#        2. Verify coordinate system definitions<br>
#        3. Review element geometry (node ordering)<br>
#        4. Test with effective modulus approach<br>
#        5. Compare isotropic vs orthotropic materials
#        </div>
       
#        <h3>Performance Issues</h3>
       
#        <h4>Slow Processing</h4>
#        <table>
#        <tr><th>Symptom</th><th>Likely Cause</th><th>Solution</th></tr>
#        <tr><td>Long connection ID phase</td><td>Complex connectivity</td><td>Reduce workers, check model complexity</td></tr>
#        <tr><td>Slow stiffness calculation</td><td>Too many workers</td><td>Reduce worker count</td></tr>
#        <tr><td>High memory usage</td><td>Large model size</td><td>Increase RAM or reduce workers</td></tr>
#        <tr><td>Process hangs</td><td>Deadlock/infinite loop</td><td>Use single worker for debugging</td></tr>
#        </table>
       
#        <h4>Memory Problems</h4>
#        <div class="formula">
#        Memory optimization strategies:
       
#        1. Reduce number of workers (less memory per worker)
#        2. Process model in batches
#        3. Clear unused data structures
#        4. Use 64-bit Python
#        5. Increase system virtual memory
#        6. Close other applications
#        </div>
       
#        <h3>Parallel Processing Issues</h3>
       
#        <h4>Worker Process Failures</h4>
#        <div class="formula">
#        Common failure modes:
       
#        1. Individual worker crashes → Isolated failure handling
#        2. Memory exhaustion → Reduce workers or batch size
#        3. Pickling errors → Check data serialization
#        4. Timeout errors → Increase process timeout
#        5. OS limits → Check system process limits
#        </div>
       
#        <h4>Inconsistent Results</h4>
#        <div class="formula">
#        Potential causes:
       
#        1. Race conditions in shared data
#        2. Non-deterministic processing order
#        3. Floating-point precision differences
#        4. Incomplete data synchronization
       
#        Solutions:
#        1. Use single worker for reproducibility testing
#        2. Check data isolation between workers
#        3. Verify deterministic algorithms
#        4. Compare multiple runs for consistency
#        </div>
       
#        <h2>Debug Mode Operation</h2>
       
#        <h3>Enabling Debug Output</h3>
#        <div class="formula">
#        Debug information includes:
       
#        1. Detailed connection analysis
#        2. Property extraction details
#        3. Coordinate transformation matrices
#        4. Intermediate calculation values
#        5. Worker process status
#        6. Memory usage statistics
#        </div>
       
#        <h3>Log Analysis</h3>
#        <div class="formula">
#        Key log entries to review:
       
#        "Connection ID": Check multi-plate grouping
#        "Property extraction": Verify material properties
#        "Stiffness calculation": Review calculated values
#        "Worker status": Monitor parallel processing
#        "Memory usage": Track resource consumption
#        </div>
       
#        <h2>Model-Specific Troubleshooting</h2>
       
#        <h3>Aerospace Models</h3>
#        <ul>
#            <li><strong>Mixed units:</strong> Ensure consistent units throughout</li>
#            <li><strong>Legacy formats:</strong> May need format conversion</li>
#            <li><strong>Large assemblies:</strong> Consider sub-model processing</li>
#            <li><strong>Complex materials:</strong> Validate composite definitions</li>
#        </ul>
       
#        <h3>Multi-Material Interfaces</h3>
#        <ul>
#            <li><strong>Metal-composite joints:</strong> Check material parameter differences</li>
#            <li><strong>Thickness mismatches:</strong>
# <li><strong>Thickness mismatches:</strong> May cause numerical issues</li>
#            <li><strong>Property discontinuities:</strong> Verify material transitions</li>
#            <li><strong>Fastener specifications:</strong> Ensure appropriate for application</li>
#        </ul>
       
#        <h3>Legacy Model Conversion</h3>
#        <ul>
#            <li><strong>Old Nastran versions:</strong> Card format differences</li>
#            <li><strong>Unit system changes:</strong> Verify consistent units</li>
#            <li><strong>Property updates:</strong> Modern material definitions</li>
#            <li><strong>Connection methods:</strong> Updated modeling practices</li>
#        </ul>
       
#        <h2>Recovery Procedures</h2>
       
#        <h3>Partial Failure Recovery</h3>
#        <div class="formula">
#        When processing partially fails:
       
#        1. Save partial results (completed elements)
#        2. Identify failed element IDs
#        3. Process failed elements separately
#        4. Merge results manually if needed
#        5. Document recovery actions
#        </div>
       
#        <h3>Data Recovery</h3>
#        <div class="formula">
#        Recovery strategies:
       
#        Backup files: Always maintain original copies
#        Checkpoint saves: Intermediate result preservation  
#        Log analysis: Reconstruct processing history
#        Manual completion: Finish remaining elements by hand
#        Alternative methods: Use different calculation approaches
#        </div>
       
#        <h2>Validation After Troubleshooting</h2>
       
#        <h3>Post-Fix Verification</h3>
#        <ol>
#            <li><strong>Rerun processing:</strong> Verify fix resolved issue</li>
#            <li><strong>Spot check results:</strong> Manual verification of samples</li>
#            <li><strong>Compare with previous:</strong> Check consistency if available</li>
#            <li><strong>Full model validation:</strong> Complete quality assurance</li>
#        </ol>
       
#        <h3>Regression Testing</h3>
#        <ul>
#            <li>Test with known good models</li>
#            <li>Verify no new issues introduced</li>
#            <li>Performance impact assessment</li>
#            <li>Documentation of changes made</li>
#        </ul>
       
#        <h2>Getting Help</h2>
       
#        <h3>Information to Collect</h3>
#        <div class="formula">
#        When seeking help, provide:
       
#        1. Complete error messages
#        2. Model characteristics (size, type, complexity)
#        3. System information (OS, RAM, CPU)
#        4. Processing parameters used
#        5. Sample problematic elements/connections
#        6. Log files (if available)
#        </div>
       
#        <h3>Support Resources</h3>
#        <ul>
#            <li><strong>Documentation:</strong> This help system and examples</li>
#            <li><strong>Log analysis:</strong> Built-in diagnostic information</li>
#            <li><strong>Test cases:</strong> Simple validation examples</li>
#            <li><strong>User community:</strong> Organization-specific support</li>
#        </ul>
       
#        <h2>Prevention Strategies</h2>
       
#        <h3>Model Quality</h3>
#        <ul>
#            <li>Use consistent modeling practices</li>
#            <li>Validate models before processing</li>
#            <li>Maintain clean, well-documented BDF files</li>
#            <li>Follow naming conventions for fastener groups</li>
#        </ul>
       
#        <h3>Processing Best Practices</h3>
#        <ul>
#            <li>Start with small test cases</li>
#            <li>Use conservative processing parameters initially</li>
#            <li>Monitor system resources during processing</li>
#            <li>Maintain backup copies of all files</li>
#        </ul>
       
#        <h3>Quality Assurance</h3>
#        <ul>
#            <li>Implement peer review processes</li>
#            <li>Use validation test suites</li>
#            <li>Document all assumptions and decisions</li>
#            <li>Maintain traceability of all changes</li>
#        </ul>
       
#        <div class="warning">
#        <strong>Important:</strong> When troubleshooting critical production models, always work with copies and maintain detailed logs of all debugging steps. Document solutions for future reference.
#        </div>
#        """
   
#     elif topic_id == "examples":
#        return base_style + """
#        <h1>Examples and Use Cases</h1>
       
#        <h2>Typical Workflow Examples</h2>
       
#        <h3>Example 1: Wing Panel Fastener Analysis</h3>
       
#        <h4>Scenario</h4>
#        <ul>
#            <li><strong>Model:</strong> Wing upper panel with 500 fasteners</li>
#            <li><strong>Materials:</strong> Carbon fiber skin, aluminum frames</li>
#            <li><strong>Fasteners:</strong> Mixed titanium bolts and Hi-Lok fasteners</li>
#            <li><strong>Goal:</strong> Update all fastener stiffness for detailed stress analysis</li>
#        </ul>
       
#        <h4>Workflow</h4>
#        <div class="formula">
#        Step 1: Load wing panel BDF file (wing_panel.bdf)
#        Step 2: Import fastener sets (wing_fastener_sets.bdf)
#        Step 3: Define fastener specifications:
#                - "NAS6204": Titanium bolt, 110 GPa
#                - "HL11V": Hi-Lok, Titanium, 110 GPa
#        Step 4: Configure processing:
#                - Workers: 6 (8-core system)
#                - Method: Normal multi-plate
#                - Material orientation: Directional
#        Step 5: Process model
#        Step 6: Export updated sets for future use
#        </div>
       
#        <h4>Expected Results</h4>
#        <ul>
#            <li>Processing time: 2-3 minutes</li>
#            <li>Success rate: >98%</li>
#            <li>Skin fasteners: Higher orientation sensitivity</li>
#            <li>Frame fasteners: More isotropic behavior</li>
#        </ul>
       
#        <h3>Example 2: Fuselage Section with Airbus Method</h3>
       
#        <h4>Scenario</h4>
#        <ul>
#            <li><strong>Model:</strong> Fuselage section with multi-plate joints</li>
#            <li><strong>Connection:</strong> Skin-frame-floor beam assemblies (4-5 plates)</li>
#            <li><strong>Requirements:</strong> Use Airbus method for complex connections</li>
#            <li><strong>Validation:</strong> Compare with existing analysis</li>
#        </ul>
       
#        <h4>Workflow</h4>
#        <div class="formula">
#        Step 1: Load fuselage model (fuselage_section.bdf)
#        Step 2: Set connection mode: "Convert to Airbus Method"
#        Step 3: Configure multi-plate method: "Normal Multi-Plate"
#        Step 4: Use default specifications for undefined fasteners
#        Step 5: Process with 4 workers
#        Step 6: Review new Airbus CBUSH elements created
#        Step 7: Export updated fastener groups including new elements
#        </div>
       
#        <h4>Results Analysis</h4>
#        <ul>
#            <li>Original connections: 120 multi-plate</li>
#            <li>New Airbus CBUSHes: 120 additional elements</li>
#            <li>Load path: Axial through Airbus, shear through regular</li>
#            <li>Validation: Consistent with previous Airbus models</li>
#        </ul>
       
#        <h2>Material-Specific Examples</h2>
       
#        <h3>Example 3: Composite Laminate Analysis</h3>
       
#        <h4>PCOMP Definition</h4>
#        <div class="formula">
#        Laminate: [0/±45/90]s IM7/8552
       
#        PCOMP   3001    0.0                             0.0
#        +       301     0.125   0.0             301     0.125   45.0
#        +       301     0.125   -45.0           301     0.125   90.0
#        +       301     0.125   90.0            301     0.125   -45.0
#        +       301     0.125   45.0            301     0.125   0.0
       
#        MAT8    301     161000. 11380.  0.32    5170.   5170.   3000.
#        </div>
       
#        <h4>CLT Results</h4>
#        <div class="formula">
#        Effective properties:
#        Ex ≈ 55000 MPa (quasi-isotropic)
#        Ey ≈ 55000 MPa
#        Gxy ≈ 20000 MPa
#        Total thickness: 1.0 mm
       
#        Fastener stiffness (6.35mm titanium):
#        K1 ≈ 3.5E7 N/mm
#        K2 ≈ K3 ≈ 2.1E6 N/mm (minimal orientation effect)
#        </div>
       
#        <h3>Example 4: Unidirectional Composite</h3>
       
#        <h4>Configuration</h4>
#        <div class="formula">
#        Laminate: [0°]8 AS4/3501-6
       
#        PCOMP   4001    0.0                             0.0
#        +       401     0.125   0.0     401     0.125   0.0
#        +       401     0.125   0.0     401     0.125   0.0
#        +       401     0.125   0.0     401     0.125   0.0
#        +       401     0.125   0.0     401     0.125   0.0
       
#        MAT8    401     142000. 10300.  0.27    7170.   7170.   4140.
#        </div>
       
#        <h4>Orientation Effects</h4>
#        <div class="formula">
#        Fastener parallel to fibers (0°):
#        Ex = 142000 MPa, Ey = 10300 MPa
#        K2 ≈ 4.2E6 N/mm (high, along fibers)
#        K3 ≈ 0.8E6 N/mm (low, across fibers)
       
#        Fastener at 45° to fibers:
#        Ex_45 ≈ 20000 MPa, Ey_45 ≈ 20000 MPa
#        K2 ≈ K3 ≈ 1.5E6 N/mm (moderate, both directions)
#        </div>
       
#        <h2>Complex Connection Examples</h2>
       
#        <h3>Example 5: Multi-Material Joint</h3>
       
#        <h4>Stack-up Configuration</h4>
#        <div class="formula">
#        Connection layers (bottom to top):
#        1. Aluminum doubler: 3.0mm, E=70000 MPa
#        2. CFRP skin: 2.0mm, [0/±45/90]s
#        3. Titanium fitting: 4.0mm, E=110000 MPa
#        4. Steel bracket: 2.5mm, E=210000 MPa
       
#        Total thickness: 11.5mm
#        Fastener: 12.7mm steel bolt
#        </div>
       
#        <h4>Analysis Approach</h4>
#        <div class="formula">
#        Method: Normal Multi-Plate
       
#        Strength analysis:
#        - Steel: 2.5 × √(210000/1000) = 36.2
#        - Titanium: 4.0 × √(110000/1000) = 132.7  ← Strongest
#        - CFRP: 2.0 × √(55000/1000) = 14.8
#        - Aluminum: 3.0 × √(70000/1000) = 25.1
       
#        → Titanium fitting selected as middle plate
#        → Average aluminum, CFRP, steel as outer equivalent
#        </div>
       
#        <h4>Results</h4>
#        <div class="formula">
#        Outer equivalent properties:
#        E_avg = (70000×3.0 + 55000×2.0 + 210000×2.5) / (3.0+2.0+2.5)
#             = 110000 MPa
#        t_avg = 7.5/3 = 2.5mm
       
#        Huth calculation:
#        Middle: Titanium, 4.0mm, 110000 MPa
#        Outer: Equivalent, 2.5mm, 110000 MPa
#        n = 2 (multi-plate factor)
#        </div>
       
#        <h3>Example 6: Repair Joint Analysis</h3>
       
#        <h4>Repair Configuration</h4>
#        <ul>
#            <li><strong>Original:</strong> 2.0mm aluminum skin</li>
#            <li><strong>Patch:</strong> 1.5mm aluminum doubler</li>
#            <li><strong>Fasteners:</strong> 4.78mm aluminum rivets</li>
#            <li><strong>Pattern:</strong> 3-row repair pattern</li>
#        </ul>
       
#        <h4>Modeling Considerations</h4>
#        <div class="formula">
#        Connection characteristics:
#        - Simple two-plate joint
#        - Identical materials (aluminum)
#        - Small fastener diameter
#        - High fastener density
       
#        Expected behavior:
#        - Symmetric stiffness (K2 = K3)
#        - Lower overall stiffness (small diameter)
#        - No orientation effects (isotropic)
#        </div>
       
#        <h2>Validation Case Studies</h2>
       
#        <h3>Case Study 1: HyperMesh Correlation</h3>
       
#        <h4>Comparison Setup</h4>
#        <div class="formula">
#        Model: Aircraft door hinge
#        Fasteners: 20× 6.35mm titanium bolts
#        Materials: Aluminum and steel components
       
#        HyperMesh results (reference):
#        Average K2 = 4.1E6 N/mm
#        Standard deviation = 3.2%
       
#        Tool results (HyperMesh method):
#        Average K2 = 4.0E6 N/mm  
#        Standard deviation = 2.8%
       
#        Correlation: 97.6% agreement
#        </div>
       
#        <h3>Case Study 2: Test Correlation</h3>
       
#        <h4>Experimental Setup</h4>
#        <ul>
#            <li><strong>Test:</strong> Lap shear specimen with single fastener</li>
#            <li><strong>Materials:</strong> 2.0mm aluminum plates</li>
#            <li><strong>Fastener:</strong> 6.35mm titanium bolt</li>
#            <li><strong>Loading:</strong> Quasi-static tension</li>
#        </ul>
       
#        <h4>Correlation Results</h4>
#        <div class="formula">
#        Test data:
#        Initial stiffness = 3.8E6 N/mm
       
#        Calculated stiffness:
#        Huth formula = 4.2E6 N/mm
#        Difference = +10.5%
       
#        Analysis:
#        - Good correlation within typical scatter
#        - Conservative prediction (higher stiffness)
#        - Accounts for initial non-linearity in test
#        </div>
       
#        <h2>Troubleshooting Examples</h2>
       
#        <h3>Example 7: Large Model Processing</h3>
       
#        <h4>Problem</h4>
#        <ul>
#            <li><strong>Model:</strong> Complete aircraft with 50,000 fasteners</li>
#            <li><strong>Issue:</strong> Memory exhaustion during processing</li>
#            <li><strong>System:</strong> 32GB RAM, 16-core workstation</li>
#        </ul>
       
#        <h4>Solution Strategy</h4>
#        <div class="formula">
#        Initial attempt: 16 workers → Memory error
       
#        Optimization steps:
#        1. Reduce workers to 8 → Still memory issues
#        2. Reduce workers to 4 → Slower but stable
#        3. Monitor memory usage → Peak 28GB
#        4. Process in sections → Break into zones
#        5. Final solution: 6 workers, zone processing
       
#        Result: Successful processing in 45 minutes
#        </div>
       
#        <h3>Example 8: Material Property Issues</h3>
       
#        <h4>Problem</h4>
#        <div class="formula">
#        Error: "Invalid MAT8 properties for material 301"
       
#        Investigation:
#        MAT8    301     150000. 9000.   0.35    4500.   4500.   3000.
       
#        Issue identified: ν₁₂ = 0.35 too high
#        Check: ν₁₂ × E₂/E₁ = 0.35 × 9000/150000 = 0.021 < 1.0 ✓
       
#        Actual issue: Tool validation too conservative
#        Solution: Accept warning, results valid
#        </div>
       
#        <h2>Best Practice Examples</h2>
       
#        <h3>Quality Assurance Workflow</h3>
#        <div class="formula">
#        Standard procedure for production models:
       
#        1. Pre-processing validation
#           - Model quality check
#           - Fastener connectivity review
#           - Material property audit
       
#        2. Processing with validation
#           - Start with 10% sample
#           - Full model processing
#           - Automated quality checks
       
#        3. Post-processing validation
#           - Statistical analysis of results
#           - Manual spot checks (5% sample)
#           - Comparison with previous analyses
#           - Engineering review and approval
#        </div>
       
#        <h3>Documentation Standards</h3>
#        <ul>
#            <li><strong>Processing log:</strong> Complete record of all operations</li>
#            <li><strong>Quality report:</strong> Statistical summary and validation</li>
#            <li><strong>Engineering notes:</strong> Assumptions and engineering judgment</li>
#            <li><strong>Approval record:</strong> Sign-off by responsible engineer</li>
#        </ul>
       
#        <div class="note">
#        <strong>Learning Path:</strong> Start with simple two-plate examples to understand the tool behavior, then progress to multi-plate and composite cases. Always validate results against known solutions or experimental data when possible.
#        </div>
#        """
   
#     else:
#         return base_style + "<h1>Topic Not Found</h1><p>The requested documentation topic could not be found.</p>"