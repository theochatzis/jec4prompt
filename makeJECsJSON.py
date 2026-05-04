import correctionlib.schemav2 as schema
import argparse
import os
import sys

# =======================================================================================
# Parser for residuals (L3Absolute & L2L3Residual) - works just with one dimension (eta)
# =======================================================================================
def parse_l2res_file(filepath, name, version=1):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().strip('{}').split()
    
    # --- DYNAMIC FORMULA EXTRACTOR ---
    try:
        corr_idx = header.index("Correction")
        
        # CMS headers always follow this pattern: 
        # { [num_bin_vars] [bin_var1...] [num_form_vars] [form_var1...] [FORMULA] Correction Name }
        num_bin_vars = int(header[0])
        num_form_vars_idx = num_bin_vars + 1
        num_form_vars = int(header[num_form_vars_idx])
        formula_start_idx = num_form_vars_idx + num_form_vars + 1
        
        # Slice the list to grab EVERY token that belongs to the formula
        formula_tokens = header[formula_start_idx : corr_idx]
        
        # Join them back together with no spaces so TFormula can read it safely
        formula_str = "".join(formula_tokens)
        
        # Fallback if the string is just a number (e.g. "1")
        if formula_str.replace('.','',1).isdigit():
            formula_str = "[0]" 
            
    except (ValueError, IndexError):
        formula_str = "[0]" 
        
    print(f"[{name}] Extracted formula: {formula_str}")

    inputs = [
        schema.Variable(name="JetEta", type="real", description="Jet Eta"),
        schema.Variable(name="JetPt", type="real", description="Jet Pt")
    ]

    eta_edges = []
    content = []
    
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"): continue
        parts = list(map(float, line.split()))
        
        eta_min, eta_max = parts[0], parts[1]
        num_params = int(parts[2])
        
        # Skip ptMin (index 3) and ptMax (index 4). Parameters start at 5.
        if num_params >= 3:
            parameters = parts[5 : 5 + (num_params - 2)]
        else:
            parameters = parts[3 : 3 + num_params]
        
        if len(eta_edges) == 0:
            eta_edges.append(eta_min)
        eta_edges.append(eta_max)
        
        formula_node = schema.Formula(
            nodetype="formula",
            expression=formula_str if formula_str != "[0]" else str(parameters[0]),
            parser="TFormula",
            variables=["JetPt"],
            parameters=parameters
        )
        content.append(formula_node)

    # Clean duplicates from eta_edges safely
    clean_edges = []
    for e in eta_edges:
        if not clean_edges or abs(e - clean_edges[-1]) > 1e-5:
            clean_edges.append(e)

    master_binning = schema.Binning(
        nodetype="binning", input="JetEta", edges=clean_edges,
        content=content, flow="clamp"
    )

    return schema.Correction(
        name=name, version=version, inputs=inputs,
        output=schema.Variable(name="Correction", type="real"), data=master_binning
    )

# =======================================================================================
# Parser for MC truth (L2Relative) - similar to the l2res but just has 2D also for Phi.
# =======================================================================================
def parse_l2_file(filepath, name="L2Relative", version=1):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().strip('{}').split()
    
    # --- Dynamic Formula Extractor ---
    try:
        correction_idx = header.index("Correction")
        formula_str = header[correction_idx - 1]
    except ValueError:
        formula_str = "[0]"
        
    print(f"[{name}] Extracted formula: {formula_str}")

    inputs = [
        schema.Variable(name="JetEta", type="real", description="Jet Eta"),
        schema.Variable(name="JetPhi", type="real", description="Jet Phi"),
        schema.Variable(name="JetPt", type="real", description="Jet Pt")
    ]
    
    eta_bins = []
    current_eta_min = None
    current_eta_max = None
    phi_content = []
    phi_edges = []
    
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"): continue
        parts = list(map(float, line.split()))
        
        eta_min, eta_max = parts[0], parts[1]
        phi_min, phi_max = parts[2], parts[3]
        num_params = int(parts[4])
        
        # CRITICAL FIX KEPT: Skip ptMin (index 5) and ptMax (index 6). 
        # Actual parameters start at 7.
        parameters = parts[7 : 7 + (num_params - 2)]
        
        formula_node = schema.Formula(
            nodetype="formula", expression=formula_str, parser="TFormula",
            variables=["JetPt"], parameters=parameters
        )

        if current_eta_min is not None and eta_min != current_eta_min:
            phi_edges.append(phi_max_prev) 
            eta_bins.append({
                "eta_min": current_eta_min, "eta_max": current_eta_max,
                "phi_edges": phi_edges, "content": phi_content
            })
            phi_content = []
            phi_edges = []
            
        if len(phi_edges) == 0:
            phi_edges.append(phi_min)
        else:
            phi_edges.append(phi_min) 
            
        phi_content.append(formula_node)
        current_eta_min = eta_min
        current_eta_max = eta_max
        phi_max_prev = phi_max

    phi_edges.append(phi_max_prev)
    eta_bins.append({
        "eta_min": current_eta_min, "eta_max": current_eta_max,
        "phi_edges": phi_edges, "content": phi_content
    })

    eta_edges = [b["eta_min"] for b in eta_bins] + [eta_bins[-1]["eta_max"]]
    eta_content = []
    
    for e_bin in eta_bins:
        clean_edges = []
        for e in e_bin["phi_edges"]:
            if not clean_edges or abs(e - clean_edges[-1]) > 1e-5:
                clean_edges.append(e)
                
        eta_content.append(schema.Binning(
            nodetype="binning", input="JetPhi", edges=clean_edges,
            content=e_bin["content"], flow="clamp"
        ))

    master_binning = schema.Binning(
        nodetype="binning", input="JetEta", edges=eta_edges,
        content=eta_content, flow="clamp"
    )

    return schema.Correction(
        name=name, version=version, inputs=inputs,
        output=schema.Variable(name="Correction", type="real"), data=master_binning
    )

# =====================================================================
# Parser for L1FastJet
# =====================================================================
def parse_l1_file(filepath, name="L1FastJet", version=1):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().strip('{}').split()
    
    # --- Dynamic Formula Extractor ---
    try:
        correction_idx = header.index("Correction")
        formula_str = header[correction_idx - 1]
    except ValueError:
        formula_str = "1"
        
    print(f"[{name}] Extracted formula: {formula_str}")

    inputs = [
        schema.Variable(name="JetEta", type="real", description="Jet Eta"),
        schema.Variable(name="JetPt", type="real", description="Jet Pt"),
        schema.Variable(name="Rho", type="real", description="Event Rho"),
        schema.Variable(name="JetA", type="real", description="Jet Area")
    ]
    
    eta_edges = []
    content = []
    
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"): continue
        parts = list(map(float, line.split()))
        
        eta_min, eta_max = parts[0], parts[1]
        num_params = int(parts[2])
        
        parameters = []
        # CRITICAL FIX KEPT: Skip the 6 limits (rhoMin, rhoMax, ptMin, ptMax, aMin, aMax).
        # Actual mathematical parameters start at index 9!
        if num_params > 6:
            parameters = parts[9 : 9 + (num_params - 6)]
        
        if len(eta_edges) == 0:
            eta_edges.append(eta_min)
        eta_edges.append(eta_max)
        
        formula_node = schema.Formula(
            nodetype="formula",
            expression=formula_str,
            parser="TFormula",
            variables=["Rho", "JetA", "JetPt"], 
            parameters=parameters
        )
        content.append(formula_node)

    clean_edges = []
    for e in eta_edges:
        if not clean_edges or abs(e - clean_edges[-1]) > 1e-5:
            clean_edges.append(e)

    master_binning = schema.Binning(
        nodetype="binning", input="JetEta", edges=clean_edges,
        content=content, flow="clamp"
    )

    return schema.Correction(
        name=name, version=version, inputs=inputs,
        output=schema.Variable(name="Correction", type="real"), data=master_binning
    )

if __name__ == "__main__":
    # Setup the argument parser
    parser = argparse.ArgumentParser(description="Convert CMS JEC .txt files into a compound correctionlib JSON.")
    
    parser.add_argument("--l1", required=True, help="Path to the L1FastJet .txt file")
    parser.add_argument("--l2", required=True, help="Path to the L2Relative .txt file")
    parser.add_argument("--l3abs", required=True, help="Path to the L3Absolute .txt file")
    parser.add_argument("--l2l3res", required=True, help="Path to the L2L3Residual .txt file")
    
    # Notice we removed the --name argument!
    parser.add_argument("--output", "-o", default="Master_JECs.json", help="Path and filename for the output JSON file")
    
    args = parser.parse_args()

    # Validate that the files actually exist
    for filepath in [args.l1, args.l2, args.l3abs, args.l2l3res]:
        if not os.path.exists(filepath):
            print(f"ERROR: Cannot find file -> {filepath}")
            sys.exit(1)

    # Automatically derive the compound name from the Residual filename
    # e.g., "Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi.txt" -> "Summer24Prompt24_V1_MC_L1L2L3Res_AK4PFPuppi"
    base_res_name = os.path.basename(args.l2l3res)
    auto_compound_name = base_res_name.replace(".txt", "").replace("L2L3Residual", "L1L2L3Res")

    print(f"--- Building JEC Payload ---")
    print(f"L1:  {args.l1}")
    print(f"L2:  {args.l2}")
    print(f"L3:  {args.l3abs}")
    print(f"Res: {args.l2l3res}")
    print(f"----------------------------")
    print(f"Derived Compound Name: {auto_compound_name}")
    
    # if MC in the base residual names make them DATA
    if 'MC' in base_res_name:
        base_res_name.replace('MC','DATA').replace(".txt", "")
    
    # Parse the files
    try: # name them based on the residual campaign version etc
        l1_corr = parse_l1_file(args.l1, base_res_name.replace(".txt", "").replace("L2L3Residual", "L1FastJet"))
        l2_corr = parse_l2_file(args.l2, base_res_name.replace(".txt", "").replace("L2L3Residual", "L2Relative"))
        l3_corr = parse_l2res_file(args.l3abs, base_res_name.replace(".txt", "").replace("L2L3Residual", "L3Absolute"))
        res_corr = parse_l2res_file(args.l2l3res, base_res_name.replace(".txt", ""))
    except Exception as e:
        print(f"\nERROR: Failed to parse one of the text files. Details: {e}")
        sys.exit(1)

    # Build the Compound Correction
    print("Gluing sequence into Compound Correction...")
    compound_jec = schema.CompoundCorrection(
        name=auto_compound_name,
        description="Full sequence applying L1, L2, L3, and Residual",
        
        inputs=[
            schema.Variable(name="JetEta", type="real", description="Jet Eta"),
            schema.Variable(name="JetPt", type="real", description="Jet Pt"),
            schema.Variable(name="JetPhi", type="real", description="Jet Phi"),
            schema.Variable(name="Rho", type="real", description="Event Rho"),
            schema.Variable(name="JetA", type="real", description="Jet Area")
        ],
        
        inputs_update=["JetPt"], 
        
        # Explicitly define how the inputs and outputs are accumulated (Multiplication)
        input_op="*",
        output_op="*",
        
        output=schema.Variable(name="Correction", type="real"),
        
        # Steps in the stack of corrections
        stack=[l1_corr.name, l2_corr.name, l3_corr.name, res_corr.name]
    )

    # Save to JSON
    print(f"Writing to output file: {args.output}")
    cset = schema.CorrectionSet(
        schema_version=2,
        description="JECs from JEC4Prompt framework.",
        
        # Individual corrections
        corrections=[l1_corr, l2_corr, l3_corr, res_corr], 
        
        # Compound corrections 
        compound_corrections=[compound_jec] 
    )

    with open(args.output, "w") as f:
        f.write(cset.model_dump_json(exclude_unset=True, indent=2))

    print("Success!")