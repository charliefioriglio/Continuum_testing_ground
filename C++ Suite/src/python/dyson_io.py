import sys
import os
import re
import numpy as np
import subprocess
import argparse
from dataclasses import dataclass
from typing import List, Sequence, Dict

# Reuse constants and logic from notebook
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
ANGULAR_LETTER_TO_L = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}

@dataclass
class ShellSpec:
    angular_momentum: int
    exponents: np.ndarray
    coefficients: np.ndarray
    is_pure: bool
    symbol: str

@dataclass
class AtomSpec:
    symbol: str
    index: int
    center_bohr: np.ndarray
    shells: List[ShellSpec]

@dataclass
class DysonInfo:
    label: str
    coefficients: np.ndarray
    transition: str
    norm: float = 1.0

@dataclass
class QChemData:
    atoms: List[AtomSpec]
    dyson_orbitals: List[DysonInfo]
    n_basis_functions: int
    pure_map: Dict[int, bool]

# --- Parsing Logic (Ported from DO_handler.ipynb) ---

def parse_float(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))

def extract_geometry(lines: Sequence[str]) -> List[AtomSpec]:
    atoms = []
    in_block = False
    idx = 0
    for line in lines:
        if "Standard Nuclear Orientation" in line:
            in_block = True
            continue
        if in_block and "----------------" in line:
            if idx > 0: return atoms # Done reading
            idx += 1 # Skip first header separator
            continue
            
        if in_block:
            tokens = line.split()
            if len(tokens) >= 5 and tokens[0].isdigit():
                symbol = tokens[1]
                # Coords in Angstroms
                x = parse_float(tokens[2]) * ANGSTROM_TO_BOHR
                y = parse_float(tokens[3]) * ANGSTROM_TO_BOHR
                z = parse_float(tokens[4]) * ANGSTROM_TO_BOHR
                atom_idx = len(atoms)
                atoms.append(AtomSpec(symbol, atom_idx, np.array([x,y,z]), []))
    return atoms

def extract_basis_block(lines: Sequence[str]) -> List[str]:
    block = []
    in_block = False
    for line in lines:
        if ("basis set" in line.lower() or "$basis" in line.lower()) and not in_block:
            # Check if it's the $rem variable or similar, we want the actual block header or $basis section
            # But $basis usually appears as a standalone line or header
            # In CuO_pVTZ.out it is '$basis'
            if "basis =" in line.lower() or "basis=" in line.lower(): continue # skip $rem
            in_block = True
            continue
        if in_block:
            if "$end" in line: break
            block.append(line)
    return block

def parse_shell(shell_lines, symbol, center_bohr, is_pure):
    header = shell_lines[0].split()
    ang_letter = header[0].upper()
    l = ANGULAR_LETTER_TO_L[ang_letter]
    n_prim = int(header[1])
    
    exps = []
    coeffs = []
    for row in shell_lines[1:]:
        parts = row.split()
        exps.append(parse_float(parts[0]))
        coeffs.append(parse_float(parts[1]))
        
    return ShellSpec(l, np.array(exps), np.array(coeffs), is_pure, symbol)

def parse_basis(basis_lines, atoms, pure_map):
    geom_iter = iter(atoms)
    current_atom = None
    cursor = 0
    
    while cursor < len(basis_lines):
        line = basis_lines[cursor].strip()
        if not line or line.startswith("$") or line.startswith("-") or line.startswith("Basis") or line.startswith("Requested"):
            cursor += 1
            continue
            
        if line == "****":
            current_atom = None
            cursor += 1
            continue
            
        tokens = line.split()
        if len(tokens) == 2 and tokens[1].isdigit() and tokens[0].isalpha():
            # Atom header
            try:
                current_atom = next(geom_iter)
            except StopIteration:
                break
            cursor += 1
            continue
            
        if current_atom is None:
            # Maybe skipping header mismatch
            cursor += 1
            continue
            
        # Shell block
        if tokens[0].upper() in ANGULAR_LETTER_TO_L:
            n_prim = int(tokens[1])
            block = [basis_lines[cursor]] + [basis_lines[cursor+i+1] for i in range(n_prim)]
            cursor += n_prim + 1
            
            l = ANGULAR_LETTER_TO_L[tokens[0].upper()]
            is_pure = pure_map.get(l, True)
            
            shell = parse_shell(block, current_atom.symbol, current_atom.center_bohr, is_pure)
            current_atom.shells.append(shell)
            continue
            
        cursor += 1
    return atoms

def total_basis_functions(atoms):
    count = 0
    for atom in atoms:
        for shell in atom.shells:
            l = shell.angular_momentum
            n_funcs = (2*l + 1) if shell.is_pure else ((l+1)*(l+2)//2)
            # n_contractions is basically 1 here per parsing logic (each P block is 1 contraction)
            # Note: SP shells are handled differently in C++ ref, but typically split in Q-Chem output?
            # Q-Chem output usually lists S and P separately unless SP used.
            # My simple parser assumes standard blocks.
            count += n_funcs
    return count

def parse_dyson(lines, n_basis):
    dysons = []
    # Simplified regex-less search for robustness
    i = 0
    current_norm = 1.0
    
    while i < len(lines):
        line = lines[i]
        
        # Capture Norm
        if "Dyson orbital norm is" in line:
            # "Left Dyson orbital norm is 0.9734"
            parts = line.split()
            try:
                current_norm = float(parts[-1])
            except:
                current_norm = 1.0
        
        if "Decomposition over AOs" in line:
            label = line.strip()
            coeffs = []
            i += 1
            while i < len(lines):
                row = lines[i].strip()
                if "Decomposition" in row or "Reference" in row or "State" in row:
                    i -= 1 # backtrack
                    break
                if row.startswith("g1p") or not row:
                    i += 1
                    continue
                
                parts = row.split()
                for p in parts:
                    if p != "*****" and not p.isalpha():
                        try:
                            coeffs.append(parse_float(p))
                        except: pass
                
                if len(coeffs) >= n_basis: break
                i += 1
            
            dysons.append(DysonInfo(label, np.array(coeffs[:n_basis]), "Unknown", current_norm))
            current_norm = 1.0 # Reset
        i += 1
    return dysons

def load_qchem(path):
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()
    
    atoms = extract_geometry(lines)
    basis_lines = extract_basis_block(lines)
    
    pure_map = {2: True, 3: True, 4: True, 5: True} # Default
    if 'purecart ="1111"' in text: pass # etc, simplistic
    # Check for Cartesian explicitly
    if 'purecart ="2222"' in text or 'purecart = 2222' in text:
        for k in pure_map: pure_map[k] = False

    # This is a critical check for G-shells
    # User had: purecart[1] == 1 -> G is pure.
    # We default pure=True.

    atoms = parse_basis(basis_lines, atoms, pure_map)
    
    # Recount basis functions to match Dyson coeffs
    n_basis = total_basis_functions(atoms)
    
    dysons = parse_dyson(lines, n_basis)
    
    return QChemData(atoms, dysons, n_basis, pure_map)

# --- Output Generation ---

def write_cpp_input(data: QChemData, dyson_indices: List[int], grid: dict, output_path: str):
    with open(output_path, "w") as f:
        atoms = data.atoms
        # Center the molecule at (0,0,0) for Cross Section calculations
        # This is crucial for valid Partial Wave Expansion (l-determination)
        coords = np.array([a.center_bohr for a in atoms])
        centroid = np.mean(coords, axis=0)
        
        # Shift atoms
        for a in atoms:
            a.center_bohr = [c - o for c, o in zip(a.center_bohr, centroid)]

        f.write(f"{len(atoms)}\n")
        for a in atoms:
            f.write(f"{a.symbol} {a.index} {a.center_bohr[0]:.6f} {a.center_bohr[1]:.6f} {a.center_bohr[2]:.6f}\n")
            
        # 2. Shells
        shells = []
        for a in atoms:
            for s in a.shells:
                shells.append((a.index, s))
                
        f.write(f"{len(shells)}\n")
        for atom_idx, s in shells:
            pure_int = 1 if s.is_pure else 0
            f.write(f"{atom_idx} {s.angular_momentum} {pure_int} {len(s.exponents)}\n")
            for e, c in zip(s.exponents, s.coefficients):
                f.write(f"{e:.6f} {c:.6f}\n")
        
        # 3. Dyson Coefficients (Supports multiple)
        f.write(f"{len(dyson_indices)}\n")
        for idx in dyson_indices:
            do = data.dyson_orbitals[idx]
            f.write(f"{len(do.coefficients)} {do.norm}\n")
            f.write(" ".join(map(str, do.coefficients)) + "\n")
            
        # 4. Grid
        f.write(f"{grid['x0']:.6f} {grid['x1']:.6f} {grid['y0']:.6f} {grid['y1']:.6f} {grid['z0']:.6f} {grid['z1']:.6f} {grid['step']:.6f}\n")

# --- Main ---

if __name__ == "__main__":
    import argparse
    import sys
    import os
    
    parser = argparse.ArgumentParser(description="Parse Q-Chem output and generate Dyson orbital grid.")
    parser.add_argument("input", help="Q-Chem output file")
    
    # Mutually exclusive selector
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dyson-index", type=int, help="Index of Single Dyson orbital (0-based)")
    group.add_argument("--dyson-pair", nargs=2, type=int, metavar=('L', 'R'), help="Indices of Left and Right Dyson orbitals")

    parser.add_argument("--grid-step", type=float, default=0.2, help="Grid step size (Bohr)")
    parser.add_argument("--padding", type=float, default=20.0, help="Padding around molecule (Bohr)")
    parser.add_argument("--cpp-bin", default="./dyson_gen", help="Path to compiled C++ binary")
    parser.add_argument("--output", default="dyson_grid.bin", help="Output binary file")
    parser.add_argument("--gen-only", action="store_true", help="Only generate the C++ input file, do not run the binary.")
    parser.add_argument("--input-out", default="cpp_input.dat", help="Name of the generated C++ input file (default: cpp_input.dat)")
    
    # Cross Section Args
    parser.add_argument("--xs", action="store_true", help="Compute Cross Sections")
    parser.add_argument("--ie", type=float, default=0.0, help="Ionization Energy (eV)")
    parser.add_argument("--lmax", type=int, default=3, help="Max Angular Momentum")
    parser.add_argument("--e-range", nargs=3, type=float, metavar=('MIN', 'MAX', 'PTS'), help="Energy range (eV): MIN MAX PTS")
    parser.add_argument("--vib-file", help="File with vibrational transitions (Energy[eV] FCF)")
    parser.add_argument("--xs-out", default="cross_section.txt", help="Output file for cross sections")

    args = parser.parse_args()
    
    print(f"Loading {args.input}...")
    data = load_qchem(args.input)
    print(f"Found {len(data.atoms)} atoms, {len(data.dyson_orbitals)} Dyson orbitals.")
    
    selected_indices = []
    if args.dyson_pair:
        i1, i2 = args.dyson_pair
        if i1 >= len(data.dyson_orbitals) or i2 >= len(data.dyson_orbitals):
            print(f"Error: Dyson indices {i1},{i2} out of range.")
            sys.exit(1)
        selected_indices = [i1, i2]
        print(f"Selected pair: {i1} (L) and {i2} (R)")
    else:
        if args.dyson_index >= len(data.dyson_orbitals):
            print(f"Error: Dyson index {args.dyson_index} out of range.")
            sys.exit(1)
        selected_indices = [args.dyson_index]
        print(f"Selected single: {args.dyson_index}")
        
    # Bounds
    coords = np.array([a.center_bohr for a in data.atoms])
    # Center the molecule at (0,0,0) for Cross Section calculations
    centroid = np.mean(coords, axis=0)
    for a in data.atoms:
        a.center_bohr = [c - o for c, o in zip(a.center_bohr, centroid)]
    
    # Recalculate bounds after centering
    coords = np.array([a.center_bohr for a in data.atoms])
    min_c = coords.min(axis=0) - args.padding
    max_c = coords.max(axis=0) + args.padding
    
    grid = {
        "x0": min_c[0], "x1": max_c[0],
        "y0": min_c[1], "y1": max_c[1],
        "z0": min_c[2], "z1": max_c[2],
        "step": args.grid_step
    }
    
    temp_inp = args.input_out
    write_cpp_input(data, selected_indices, grid, temp_inp)
    
    if args.xs:
        if not args.e_range:
            print("Error: --xs requires --e-range MIN MAX PTS")
            sys.exit(1)
        
        # Read Vibrational Data if provided
        vib_states = [] # (E_bind, FCF)
        if args.vib_file:
            print(f"Reading vibrational data from {args.vib_file}...")
            try:
                with open(args.vib_file, "r") as f:
                    for line in f:
                        if line.strip() and not line.startswith("#"):
                            parts = line.split()
                            if len(parts) >= 2:
                                vib_states.append((float(parts[0]), float(parts[1])))
            except Exception as e:
                print(f"Error reading vib file: {e}")
                sys.exit(1)
            print(f"Found {len(vib_states)} vibrational states.")
            # Use the first state as the base IE for C++ calculation to generate the curve
            # ACTUALLY: The C++ code takes IE as input. 
            # Strategy:
            # 1. We need sigma_el(eKE).
            # 2. C++ computes sigma_el(E_ph - IE).
            # 3. We should run C++ with IE=0.0 to get sigma_el(E_ph) = sigma_el(eKE).
            # 4. Then in python we map: sigma_v(E) = sigma_el(E - E_bind_v) * FCF^2.
            
            # So, run C++ with --ie 0.0 and range [0, Max_E_needed].
            # Max E_ph user asks for is e_range[1]. 
            # Max eKE needed is e_range[1] - min(E_bind).
            # Min eKE needed is 0.
            
            min_bind = min(v[0] for v in vib_states)
            max_e_ph = args.e_range[1]
            max_eKE = max_e_ph # Assuming IE=0
            
            # Override IE passed to C++ to be 0 for purely electronic calc
            # But wait, C++ might have 1/k checks or threshold behavior.
            # Best to let C++ just compute sigma(E) for a range of Energies which are interpreted as eKE.
            # If we pass IE=0 to C++, then E_ph = eKE.
            
            print("Computing base electronic cross section (IE=0 effectively)...")
            ie_for_cpp = 0.0
            min_eKE = 0.01 # Avoid 0
            max_eKE = args.e_range[1] - min_bind + 0.5 # Add buffer
            n_pts_cpp = 200 # High res for interpolation
            
            with open(temp_inp, "a") as f:
                f.write(f"\n{ie_for_cpp} {args.lmax} {n_pts_cpp} {min_eKE} {max_eKE}\n")
                
        else:
            # Standard single channel calculation
            if args.ie <= 0:
                print("Warning: IE is 0 or negative. Ensure this is correct.")
            with open(temp_inp, "a") as f:
                # User specified e-range as eKE (if IE provided)
                # C++ expects Photon Energy range.
                # E_ph = eKE + IE
                e_min_ph = args.e_range[0] + args.ie
                e_max_ph = args.e_range[1] + args.ie
                f.write(f"\n{args.ie} {args.lmax} {int(args.e_range[2])} {e_min_ph} {e_max_ph}\n")
    
    print(f"Input file written to {temp_inp}")
    
    if args.gen_only:
        print("Generation only requested. Exiting.")
        sys.exit(0)

    print(f"Running {args.cpp_bin}...")
    if not os.path.exists(args.cpp_bin):
        print(f"Error: C++ binary {args.cpp_bin} not found. Did you compile?")
        print("Run: make dyson_gen")
        sys.exit(1)
        
    subprocess.run([args.cpp_bin, temp_inp, args.output], check=True)
    
    # Post-processing for Relative Cross Sections
    if args.xs and args.vib_file:
         print("Processing relative cross sections...")
         # Read C++ output: eKE (since IE=0) vs Sigma
         # cross_section.txt format: Energy(eV) CrossSection(au)
         try:
             raw_data = np.loadtxt("cross_section.txt") # This is E_ph(=eKE) vs Sigma
             eKE_grid = raw_data[:, 0]
             sigma_el_grid = raw_data[:, 1]
         except Exception as e:
             print(f"Error reading cross_section.txt: {e}")
             sys.exit(1)
             
         from scipy.interpolate import interp1d
         # Create interpolator for electronic cross section
         # Fill value 0 for eKE < 0 or > max
         sigma_func = interp1d(eKE_grid, sigma_el_grid, kind='cubic', bounds_error=False, fill_value=0.0)
         
         # User requested E_ph grid
         user_e_min, user_e_max, user_n_pts = args.e_range
         E_ph_list = np.linspace(user_e_min, user_e_max, int(user_n_pts))
         
         total_sigma = np.zeros_like(E_ph_list)
         channel_sigmas = []
         
         for (E_bind, fcf) in vib_states:
             # Calculate eKE for this channel at each photon energy
             # eKE = E_ph - E_bind
             ke_vals = E_ph_list - E_bind
             # sigma_v = sigma_el(ke) * |FCF|^2
             sig_v = sigma_func(ke_vals) * (fcf**2)
             # Ensure no negative cross sections from spline overshoot
             sig_v = np.maximum(sig_v, 0.0)
             
             total_sigma += sig_v
             channel_sigmas.append(sig_v)
             
         # Write output
         out_file = "xs_relative.csv"
         print(f"Writing relative cross sections to {out_file}...")
         header = "E_photon(eV),Total_XS"
         for i in range(len(vib_states)):
             header += f",Ch{i}_Ebind={vib_states[i][0]:.3f}"
             
         out_data = np.column_stack([E_ph_list, total_sigma] + channel_sigmas)
         np.savetxt(out_file, out_data, delimiter=",", header=header, comments="")
         
    print("Done.")
