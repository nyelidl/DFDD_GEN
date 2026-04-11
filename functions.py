#@title **LB-PaCS-MD Backend (`functions.py`)** ⚙️
#@markdown Core functions for LB-PaCS-MD: PDB cleanup, validation, protonation, ligand setup, tleap, minimization, and output packaging.

import os
import re
import math
import shutil
import subprocess
import urllib.request
import zipfile
import json
import warnings
from pathlib import Path
from datetime import datetime
from typing import Optional

import numpy as np

# ============================================================
# CONSTANTS
# ============================================================

SILENT_DEFAULTS = {
    "min_max_iterations":    5000,
    "min_tolerance":         10.0,
    "heat_stages":           3,
    "heat_steps_per_stage":  5000,
    "nvt_steps":             5000,
    "timestep_fs":           2.0,
    "ntc":                   2,
    "ntf":                   2,
    "cut":                   10.0,
    "nonbonded_method":      "PME",
    "ewald_tolerance":       0.0005,
    "ntp":                   1,
    "barostat":              2,
    "pres0":                 1.0,
    "ntt":                   3,
    "gamma_ln":              2.0,
    "iwrap":                 1,
    "ntwx":                  500,
    "box_type":              "octahedral",
    "box_padding":           15.0,
    "ions":                  "Na+ 0, Cl- 0",
    "pacs_mode":             "-b",
    "bsv_sphere_radius":     20.0,
    "bsv_subdivisions":      3,
    "bsv_hemisphere":        True,
    "min_placement_dist":    30.0,
    "default_buffer":        35.0,
}

# HETATM classification sets
AMBER_MODIFIED_RESIDUES = {
    "PTR", "SEP", "TPO", "HYP", "MLY", "MLZ",
    "M3L", "HYL", "ALY", "ACE", "NME", "CSO",
    "CSD", "OCS", "KCX", "LYZ", "MSE",
}

KNOWN_COFACTORS = {
    "NAD", "NAI", "FAD", "FMN", "ATP", "ADP",
    "AMP", "COA", "PLP",
}

CRYSTAL_ADDITIVES_REMOVE = {
    "GOL", "PEG", "MPD", "EDO", "SO4", "PO4",
    "ACT", "CIT", "BME", "DTT", "TCE", "DMS",
    "ACN", "MES", "FMT",
}

DETERGENTS_REMOVE = {
    "DDM", "OG", "LDAO", "C8E4", "NG", "DM", "CHAPS", "SDS",
}

IONS_KEEP = {"ZN", "MG", "CA", "MN", "FE", "CU", "NI", "CO"}

IONS_REMOVE = {"NA", "CL", "K", "BR", "IOD"}

CONVERT_AUTO = {"MSE": "MET"}

NUCLEIC_RESNAMES = {
    "DA", "DT", "DG", "DC", "A", "U", "G", "C",
    "RA", "RU", "RG", "RC", "ADE", "CYT", "GUA",
    "THY", "URA",
}

LIPID_RESNAMES = {
    "POPC", "POPE", "POPG", "CHOL", "DPPC",
    "DPPE", "DPPG", "DMPC", "DLPC", "PLPC",
    "CHL1", "POPS",
}

METAL_CLUSTER_RESNAMES = {
    "HEM", "HEA", "HEB", "HEC", "FES", "CUA",
    "CUB", "NFS", "F3S", "SF4", "MOS",
}

WATER_RESNAMES = {"HOH", "WAT", "TIP", "SOL", "H2O", "TP3"}

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75,
    "BR": 1.85, "I":  1.98, "ZN": 1.39, "MG": 1.73,
    "CA": 1.49, "MN": 1.61, "FE": 1.56, "CU": 1.40,
    "NI": 1.63, "CO": 1.26,
}


class HardStopError(Exception):
    """Raised when the system cannot be processed."""
    pass


class WorkflowWarning(UserWarning):
    """Non-fatal workflow warning."""
    pass


# ============================================================
# SECTION A — PDB LOADING & CLEANING
# ============================================================

def load_pdb(source: str) -> dict:
    """
    Load a PDB from RCSB ID or local file path.

    Parameters
    ----------
    source : str
        4-character PDB ID (e.g. "1HSG") or path to .pdb file.

    Returns
    -------
    dict with keys:
        pdb_path   : str   — path to saved .pdb file
        pdb_id     : str   — PDB ID if downloaded, else filename stem
        lines      : list  — raw lines of PDB file
        warnings   : list
        errors     : list
    """
    result = {"pdb_path": None, "pdb_id": None,
              "lines": [], "warnings": [], "errors": []}

    source = source.strip()

    # Determine if RCSB ID or file path
    if len(source) == 4 and source.isalnum() and not os.path.exists(source):
        pdb_id   = source.upper()
        pdb_path = f"{pdb_id}.pdb"
        url      = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            urllib.request.urlretrieve(url, pdb_path)
        except Exception as e:
            result["errors"].append(f"Failed to download {pdb_id}: {e}")
            return result
        result["pdb_id"] = pdb_id
    else:
        if not os.path.exists(source):
            result["errors"].append(f"File not found: {source}")
            return result
        pdb_path = source
        result["pdb_id"] = Path(source).stem

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    result["pdb_path"] = pdb_path
    result["lines"]    = lines
    return result


def strip_conect_altloc(pdb_path: str, output_path: str = None) -> dict:
    """
    Remove CONECT records and resolve alternate locations (keep highest occupancy).

    Parameters
    ----------
    pdb_path    : str — input PDB path
    output_path : str — output path (default: pdb_path with _clean suffix)

    Returns
    -------
    dict with keys: pdb_path, warnings, errors
    """
    result = {"pdb_path": None, "warnings": [], "errors": []}

    if output_path is None:
        stem        = Path(pdb_path).stem
        output_path = f"{stem}_clean.pdb"

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Track best occupancy per (chain, resid, icode, atomname)
    atom_best   = {}   # key → (occupancy, line)
    kept_lines  = []
    altloc_seen = False

    for line in lines:
        rec = line[:6].strip()

        # Skip CONECT entirely
        if rec == "CONECT":
            continue

        if rec in ("ATOM", "HETATM"):
            altloc   = line[16]
            chain    = line[21]
            resseq   = line[22:26].strip()
            icode    = line[26]
            atomname = line[12:16].strip()
            key      = (chain, resseq, icode, atomname)

            if altloc not in (" ", "A", ""):
                altloc_seen = True
                try:
                    occ = float(line[54:60].strip() or "1.0")
                except ValueError:
                    occ = 1.0

                if key not in atom_best or occ > atom_best[key][0]:
                    # Blank out altloc field
                    clean = line[:16] + " " + line[17:]
                    atom_best[key] = (occ, clean)
                continue

            elif altloc == "A":
                # Explicit altloc A — treat as primary
                try:
                    occ = float(line[54:60].strip() or "1.0")
                except ValueError:
                    occ = 1.0
                clean = line[:16] + " " + line[17:]
                if key not in atom_best or occ > atom_best[key][0]:
                    atom_best[key] = (occ, clean)
                continue

        kept_lines.append(line)

    # Inject resolved altloc atoms before END/TER
    resolved = [v[1] for v in atom_best.values()]

    # Rebuild: insert resolved atoms after last ATOM/HETATM in kept_lines
    final_lines = []
    atom_block  = []
    tail_lines  = []
    in_tail     = False

    for line in kept_lines:
        rec = line[:6].strip()
        if rec in ("END", "TER", "MASTER"):
            in_tail = True
        if in_tail:
            tail_lines.append(line)
        else:
            atom_block.append(line)

    final_lines = atom_block + resolved + tail_lines

    with open(output_path, "w") as f:
        f.writelines(final_lines)

    if altloc_seen:
        result["warnings"].append(
            "Alternate locations detected — kept highest occupancy only."
        )

    result["pdb_path"] = output_path
    return result


def check_biological_assembly(pdb_path: str) -> dict:
    """
    Parse REMARK 350 to check if the biological assembly
    differs from the asymmetric unit.

    Returns
    -------
    dict with keys:
        is_multimer  : bool
        n_chains     : int
        remark_350   : str  — raw REMARK 350 text
        warnings     : list
        errors       : list
    """
    result = {
        "is_multimer": False,
        "n_chains":    1,
        "remark_350":  "",
        "warnings":    [],
        "errors":      [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    remark_lines = [l for l in lines if l.startswith("REMARK 350")]
    result["remark_350"] = "".join(remark_lines)

    # Look for "APPLY THE FOLLOWING TO CHAINS:"
    chains_in_assembly = []
    for line in remark_lines:
        if "APPLY THE FOLLOWING TO CHAINS:" in line:
            chain_part = line.split("CHAINS:")[-1].strip().rstrip(",")
            chains_in_assembly += [c.strip() for c in chain_part.split(",")]

    if len(chains_in_assembly) > 1:
        result["is_multimer"] = True
        result["n_chains"]    = len(chains_in_assembly)
        result["warnings"].append(
            f"REMARK 350: Biological assembly has {len(chains_in_assembly)} chains "
            f"({', '.join(chains_in_assembly)}). "
            "If only one chain is the target, use filter_chains()."
        )

    return result


def filter_chains(pdb_path: str, chains: list, output_path: str = None) -> dict:
    """
    Keep only specified chains from PDB.

    Parameters
    ----------
    pdb_path    : str  — input PDB
    chains      : list — e.g. ["A"] or ["A", "B"]
    output_path : str  — default: stem_chainX.pdb

    Returns
    -------
    dict with keys: pdb_path, kept_chains, warnings, errors
    """
    result = {"pdb_path": None, "kept_chains": chains,
              "warnings": [], "errors": []}

    if output_path is None:
        chain_tag   = "".join(chains)
        output_path = f"{Path(pdb_path).stem}_chain{chain_tag}.pdb"

    chains_set = set(c.upper() for c in chains)

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    kept = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM", "TER"):
            if len(line) > 21 and line[21].upper() in chains_set:
                kept.append(line)
        elif rec in ("END", "MASTER", "REMARK", "HEADER",
                     "TITLE", "SEQRES", "CRYST1"):
            kept.append(line)

    if not kept:
        result["errors"].append(
            f"No atoms found for chains {chains}. Check chain IDs in PDB."
        )
        return result

    with open(output_path, "w") as f:
        f.writelines(kept)
        if not kept[-1].startswith("END"):
            f.write("END\n")

    result["pdb_path"] = output_path
    return result


def renumber_residues(pdb_path: str, output_path: str = None,
                      start: int = 1) -> dict:
    """
    Renumber residues sequentially starting from `start`.
    Handles multiple chains independently.

    Returns
    -------
    dict with keys: pdb_path, resid_map, warnings, errors
    resid_map: {(chain, old_resid): new_resid}
    """
    result = {"pdb_path": None, "resid_map": {},
              "warnings": [], "errors": []}

    if output_path is None:
        output_path = f"{Path(pdb_path).stem}_renumbered.pdb"

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Build mapping per chain
    chain_residues = {}
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            chain  = line[21]
            resseq = line[22:26].strip()
            icode  = line[26]
            key    = (chain, resseq, icode)
            if chain not in chain_residues:
                chain_residues[chain] = []
            if key not in chain_residues[chain]:
                chain_residues[chain].append(key)

    # Build sequential map
    resid_map = {}
    for chain, residues in chain_residues.items():
        for new_num, (ch, old_seq, icode) in enumerate(residues, start=start):
            resid_map[(ch, old_seq, icode)] = new_num

    # Apply map
    new_lines = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            chain  = line[21]
            resseq = line[22:26].strip()
            icode  = line[26]
            key    = (chain, resseq, icode)
            if key in resid_map:
                new_num  = resid_map[key]
                new_line = line[:22] + f"{new_num:4d}" + line[26:]
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(new_lines)

    # Simplify resid_map for return: {(chain, old_resid): new_resid}
    simple_map = {(ch, old): new
                  for (ch, old, ic), new in resid_map.items()}
    result["pdb_path"]  = output_path
    result["resid_map"] = simple_map
    return result


def convert_mse_to_met(pdb_path: str, output_path: str = None) -> dict:
    """
    Convert selenomethionine (MSE) to methionine (MET).
    Renames residue and changes SE → SD atom name.

    Returns
    -------
    dict with keys: pdb_path, n_converted, warnings, errors
    """
    result = {"pdb_path": None, "n_converted": 0,
              "warnings": [], "errors": []}

    if output_path is None:
        output_path = f"{Path(pdb_path).stem}_mse2met.pdb"

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    count = 0
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            resname = line[17:20].strip()
            if resname == "MSE":
                atomname = line[12:16].strip()
                new_atom = " SD " if atomname == "SE" else line[12:16]
                new_line = (line[:6]
                            + line[6:12]
                            + f"{new_atom:4s}"
                            + line[16:17]
                            + "MET"
                            + line[20:])
                new_lines.append(new_line)
                count += 1
                continue
        new_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(new_lines)

    if count > 0:
        result["warnings"].append(
            f"Converted {count} MSE atoms to MET (SE → SD)."
        )
    result["pdb_path"]    = output_path
    result["n_converted"] = count
    return result


# ============================================================
# SECTION B — HARD STOP CHECKS
# ============================================================

def check_hard_stops(pdb_path: str) -> dict:
    """
    Scan PDB for system types that cannot be handled by this workflow.
    Raises HardStopError immediately if found.

    Checks:
        - Nucleic acids (DNA/RNA)
        - Membrane lipids
        - Heme / Fe-S clusters / metal clusters

    Returns
    -------
    dict with keys: passed, warnings, errors
    """
    result = {"passed": False, "warnings": [], "errors": []}

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    found_nucleic = set()
    found_lipids  = set()
    found_metals  = set()

    for line in lines:
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = line[17:20].strip().upper()

        if resname in NUCLEIC_RESNAMES:
            found_nucleic.add(resname)
        if resname in LIPID_RESNAMES:
            found_lipids.add(resname)
        if resname in METAL_CLUSTER_RESNAMES:
            found_metals.add(resname)

    messages = []

    if found_nucleic:
        messages.append(
            f"NUCLEIC ACIDS detected: {found_nucleic}\n"
            "  → This workflow supports protein-ligand systems only.\n"
            "  → For DNA/RNA systems, use OL15/OL3 force fields manually."
        )

    if found_lipids:
        messages.append(
            f"MEMBRANE LIPIDS detected: {found_lipids}\n"
            "  → Membrane protein systems are not supported.\n"
            "  → Recommended: use CHARMM-GUI Membrane Builder."
        )

    if found_metals:
        messages.append(
            f"HEME / METAL CLUSTERS detected: {found_metals}\n"
            "  → Complex metal centers require specialized parameterization.\n"
            "  → Recommended: use MCPB.py for metal cluster setup."
        )

    if messages:
        full_msg = "\n\n".join(["❌ HARD STOP — System cannot be processed:"]
                               + messages)
        raise HardStopError(full_msg)

    result["passed"] = True
    return result


# ============================================================
# SECTION C — STRUCTURE VALIDATION + HETATM CLASSIFICATION
# ============================================================

def check_missing_residues(pdb_path: str) -> dict:
    """
    Detect missing residues from REMARK 465.
    Raises HardStopError if missing residues found in the protein chain.

    Returns
    -------
    dict with keys: missing_residues, passed, warnings, errors
    missing_residues: list of dicts {chain, resname, resid}
    """
    result = {
        "missing_residues": [],
        "passed":           False,
        "warnings":         [],
        "errors":           [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    missing = []
    in_465  = False

    for line in lines:
        if "REMARK 465" in line:
            in_465 = True
            # Lines look like: REMARK 465   MET A     1
            parts = line[10:].strip().split()
            if len(parts) >= 3:
                resname = parts[0]
                chain   = parts[1] if len(parts[1]) == 1 else "?"
                try:
                    resid = int(parts[2])
                    missing.append({
                        "chain":   chain,
                        "resname": resname,
                        "resid":   resid,
                    })
                except (ValueError, IndexError):
                    pass
        elif in_465 and not line.startswith("REMARK"):
            break

    if missing:
        resid_list = ", ".join(
            f"{m['chain']}:{m['resname']}{m['resid']}" for m in missing[:10]
        )
        if len(missing) > 10:
            resid_list += f" ... (+{len(missing)-10} more)"

        raise HardStopError(
            f"❌ HARD STOP — Missing residues detected ({len(missing)} total):\n"
            f"  {resid_list}\n\n"
            "  Missing residues must be modeled before MD simulation.\n"
            "  Recommended tools:\n"
            "    • Modeller: https://salilab.org/modeller/\n"
            "    • RoseTTAFold: https://robetta.bakerlab.org/\n"
            "    • AlphaFold2 + Phenix for loop building"
        )

    result["passed"]           = True
    result["missing_residues"] = missing
    return result


def check_chain_breaks(pdb_path: str,
                       max_ca_dist: float = 4.0) -> dict:
    """
    Detect chain breaks by checking Cα–Cα distances between consecutive residues.
    Expected Cα–Cα distance for peptide bond ≈ 3.8 Å.

    Parameters
    ----------
    max_ca_dist : float — threshold above which a break is flagged (default 4.0 Å)

    Returns
    -------
    dict with keys: breaks, passed, warnings, errors
    breaks: list of dicts {chain, resid_before, resid_after, distance}
    """
    result = {
        "breaks":   [],
        "passed":   False,
        "warnings": [],
        "errors":   [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Collect CA atoms per chain
    ca_per_chain = {}
    for line in lines:
        rec = line[:6].strip()
        if rec != "ATOM":
            continue
        atomname = line[12:16].strip()
        if atomname != "CA":
            continue
        chain  = line[21]
        resseq = int(line[22:26].strip())
        x      = float(line[30:38].strip())
        y      = float(line[38:46].strip())
        z      = float(line[46:54].strip())

        if chain not in ca_per_chain:
            ca_per_chain[chain] = []
        ca_per_chain[chain].append((resseq, np.array([x, y, z])))

    breaks = []
    for chain, ca_list in ca_per_chain.items():
        ca_sorted = sorted(ca_list, key=lambda x: x[0])
        for i in range(len(ca_sorted) - 1):
            res_a, coord_a = ca_sorted[i]
            res_b, coord_b = ca_sorted[i + 1]

            # Only flag consecutive residues (gap of 1)
            if res_b - res_a != 1:
                continue

            dist = float(np.linalg.norm(coord_b - coord_a))
            if dist > max_ca_dist:
                breaks.append({
                    "chain":        chain,
                    "resid_before": res_a,
                    "resid_after":  res_b,
                    "distance":     round(dist, 2),
                })

    if breaks:
        break_list = ", ".join(
            f"{b['chain']}:{b['resid_before']}-{b['resid_after']} "
            f"({b['distance']:.1f} Å)"
            for b in breaks[:5]
        )
        raise HardStopError(
            f"❌ HARD STOP — Chain breaks detected ({len(breaks)} total):\n"
            f"  {break_list}\n\n"
            "  Chain breaks must be resolved before simulation.\n"
            "  Possible causes: missing loops, crystal packing artifacts.\n"
            "  Recommended: use Modeller or AlphaFold2 to fill gaps."
        )

    result["passed"] = True
    result["breaks"] = breaks
    return result


def classify_hetatm(pdb_path: str) -> dict:
    """
    Classify all HETATM records into categories.

    Categories:
        amber_modified    : modified residues handled by ff19SB
        cofactors         : known cofactors (NAD, FAD, etc.)
        crystal_additives : remove automatically
        detergents        : remove automatically
        ions_keep         : metal ions (ZN, MG, etc.) — keep
        ions_remove       : buffer ions (NA, CL, etc.) — remove (tleap re-adds)
        waters            : crystallographic waters
        convert           : MSE → MET auto-conversion
        unknown           : ask user
        cocrystal_ligands : potential drug/ligand molecules

    Returns
    -------
    dict with full classification per resname
    """
    result = {
        "amber_modified":    [],
        "cofactors":         [],
        "crystal_additives": [],
        "detergents":        [],
        "ions_keep":         [],
        "ions_remove":       [],
        "waters":            [],
        "convert":           [],
        "unknown":           [],
        "cocrystal_ligands": [],
        "warnings":          [],
        "errors":            [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    seen = {}
    for line in lines:
        rec = line[:6].strip()
        if rec != "HETATM":
            continue
        resname = line[17:20].strip().upper()
        chain   = line[21]
        resseq  = line[22:26].strip()
        key     = (resname, chain, resseq)
        if key in seen:
            continue
        seen[key] = True

        entry = {"resname": resname, "chain": chain, "resid": resseq}

        if resname in WATER_RESNAMES:
            result["waters"].append(entry)
        elif resname in AMBER_MODIFIED_RESIDUES:
            result["amber_modified"].append(entry)
        elif resname in KNOWN_COFACTORS:
            result["cofactors"].append(entry)
        elif resname in CRYSTAL_ADDITIVES_REMOVE:
            result["crystal_additives"].append(entry)
        elif resname in DETERGENTS_REMOVE:
            result["detergents"].append(entry)
        elif resname in IONS_KEEP:
            result["ions_keep"].append(entry)
        elif resname in IONS_REMOVE:
            result["ions_remove"].append(entry)
        elif resname in CONVERT_AUTO:
            result["convert"].append(entry)
        else:
            # Heuristic: likely a drug/ligand if MW > small ion
            # We'll flag as cocrystal_ligand if it has many atoms
            result["cocrystal_ligands"].append(entry)

    return result


def check_nonstandard_residues(pdb_path: str) -> dict:
    """
    List all non-standard residues in ATOM records
    (residues not in the standard 20 amino acids).

    Returns
    -------
    dict with keys: nonstandard, warnings, errors
    """
    STANDARD_AA = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
        "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
        "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        # Protonation variants
        "HID", "HIE", "HIP", "HIN",
        "CYX", "CYM",
        "ASH", "GLH", "LYN",
    }

    result = {"nonstandard": [], "warnings": [], "errors": []}

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    seen = set()
    for line in lines:
        rec = line[:6].strip()
        if rec != "ATOM":
            continue
        resname = line[17:20].strip().upper()
        chain   = line[21]
        resseq  = line[22:26].strip()
        key     = (resname, chain, resseq)

        if resname not in STANDARD_AA and key not in seen:
            seen.add(key)
            result["nonstandard"].append({
                "resname": resname,
                "chain":   chain,
                "resid":   resseq,
            })

    if result["nonstandard"]:
        names = list({r["resname"] for r in result["nonstandard"]})
        result["warnings"].append(
            f"Non-standard residues in ATOM records: {names}. "
            "Verify these are handled correctly by the force field."
        )

    return result


# ============================================================
# SECTION D — BINDING SITE DETECTION
# ============================================================

def detect_cocrystal_ligand(pdb_path: str,
                             hetatm_classification: dict = None) -> dict:
    """
    Detect the co-crystal ligand from HETATM records.
    Uses the cocrystal_ligands list from classify_hetatm().
    If multiple candidates, returns all and flags for user selection.

    Returns
    -------
    dict with keys:
        candidates       : list of {resname, chain, resid, n_atoms}
        selected         : dict — auto-selected if only one candidate
        warnings         : list
        errors           : list
    """
    result = {
        "candidates": [],
        "selected":   None,
        "warnings":   [],
        "errors":     [],
    }

    if hetatm_classification is None:
        hetatm_classification = classify_hetatm(pdb_path)

    candidates = hetatm_classification.get("cocrystal_ligands", [])

    if not candidates:
        result["errors"].append(
            "No co-crystal ligand found in HETATM records.\n"
            "Ensure the PDB file contains the ligand as HETATM records."
        )
        return result

    # Count atoms per candidate
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    atom_counts = {}
    for line in lines:
        rec = line[:6].strip()
        if rec != "HETATM":
            continue
        resname = line[17:20].strip().upper()
        chain   = line[21]
        resseq  = line[22:26].strip()
        key     = (resname, chain, resseq)
        atom_counts[key] = atom_counts.get(key, 0) + 1

    enriched = []
    for cand in candidates:
        key = (cand["resname"], cand["chain"], cand["resid"])
        enriched.append({
            **cand,
            "n_atoms": atom_counts.get(key, 0),
        })

    # Filter out tiny molecules (< 5 heavy atoms → likely ions or fragments)
    drug_like = [c for c in enriched if c["n_atoms"] >= 5]

    if not drug_like:
        drug_like = enriched  # fallback: keep all

    result["candidates"] = drug_like

    if len(drug_like) == 1:
        result["selected"] = drug_like[0]
    else:
        result["warnings"].append(
            f"Multiple ligand candidates found: "
            f"{[c['resname'] for c in drug_like]}. "
            "User must select the target ligand."
        )

    return result


def detect_binding_site(pdb_path: str,
                        ligand_resname: str,
                        ligand_chain: str = None,
                        ligand_resid:  str = None,
                        cutoff: float = 4.0) -> dict:
    """
    Detect binding site residues within `cutoff` Å of the co-crystal ligand.
    Only considers protein ATOM records (standard amino acids).

    Parameters
    ----------
    pdb_path       : str   — path to PDB containing both protein and ligand
    ligand_resname : str   — residue name of ligand (e.g. "LIG")
    ligand_chain   : str   — chain of ligand (optional filter)
    ligand_resid   : str   — resid of ligand (optional filter)
    cutoff         : float — distance cutoff in Å (default 4.0)

    Returns
    -------
    dict with keys:
        binding_site_resids : list of {chain, resid, resname, tleap_resid}
        binding_site_center : np.ndarray (3,) — centroid of Cα atoms
        n_residues          : int
        warnings            : list
        errors              : list
    """
    result = {
        "binding_site_resids": [],
        "binding_site_center": None,
        "n_residues":          0,
        "warnings":            [],
        "errors":              [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Collect ligand heavy atom coordinates
    lig_coords = []
    for line in lines:
        rec = line[:6].strip()
        if rec != "HETATM":
            continue
        resname = line[17:20].strip().upper()
        chain   = line[21]
        resseq  = line[22:26].strip()

        if resname != ligand_resname.upper():
            continue
        if ligand_chain and chain != ligand_chain:
            continue
        if ligand_resid and resseq != str(ligand_resid):
            continue

        # Skip hydrogens
        element = line[76:78].strip().upper() if len(line) > 76 else ""
        atomname = line[12:16].strip()
        if element == "H" or atomname.startswith("H"):
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            lig_coords.append(np.array([x, y, z]))
        except ValueError:
            continue

    if not lig_coords:
        result["errors"].append(
            f"No heavy atoms found for ligand {ligand_resname}. "
            "Check resname, chain, and resid."
        )
        return result

    lig_coords_arr = np.array(lig_coords)

    # Collect protein residue atoms
    prot_atoms = {}   # (chain, resseq, resname) → list of coords
    ca_coords  = {}   # (chain, resseq) → CA coord

    for line in lines:
        rec = line[:6].strip()
        if rec != "ATOM":
            continue
        chain    = line[21]
        resseq   = line[22:26].strip()
        resname  = line[17:20].strip().upper()
        atomname = line[12:16].strip()

        element  = line[76:78].strip().upper() if len(line) > 76 else ""
        if element == "H" or atomname.startswith("H"):
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue

        key = (chain, resseq, resname)
        if key not in prot_atoms:
            prot_atoms[key] = []
        prot_atoms[key].append(np.array([x, y, z]))

        if atomname == "CA":
            ca_coords[(chain, resseq)] = np.array([x, y, z])

    # Find residues within cutoff of any ligand atom
    binding_residues = []
    for (chain, resseq, resname), coords in prot_atoms.items():
        in_site = False
        for pcoord in coords:
            dists = np.linalg.norm(lig_coords_arr - pcoord, axis=1)
            if np.min(dists) <= cutoff:
                in_site = True
                break
        if in_site:
            binding_residues.append({
                "chain":      chain,
                "resid":      resseq,
                "resname":    resname,
                "tleap_resid": None,  # filled by build_tleap_resid_map()
            })

    if not binding_residues:
        result["warnings"].append(
            f"No protein residues found within {cutoff} Å of {ligand_resname}. "
            "Try increasing the cutoff."
        )
        return result

    # Calculate binding site center from Cα atoms
    ca_list = []
    for res in binding_residues:
        key = (res["chain"], res["resid"])
        if key in ca_coords:
            ca_list.append(ca_coords[key])

    if ca_list:
        center = np.mean(ca_list, axis=0)
    else:
        # Fallback: centroid of all binding site atom coords
        all_coords = []
        for res in binding_residues:
            key = (res["chain"], res["resid"], res["resname"])
            all_coords.extend(prot_atoms.get(key, []))
        center = np.mean(all_coords, axis=0) if all_coords else np.zeros(3)

    result["binding_site_resids"] = binding_residues
    result["binding_site_center"] = center
    result["n_residues"]          = len(binding_residues)
    return result


def detect_binding_site_multichain(pdb_path: str,
                                   ligand_resname: str,
                                   cutoff: float = 4.0) -> dict:
    """
    Wrapper for detect_binding_site that handles multi-chain proteins.
    Searches across all chains.

    Returns same dict as detect_binding_site().
    """
    return detect_binding_site(
        pdb_path       = pdb_path,
        ligand_resname = ligand_resname,
        ligand_chain   = None,
        ligand_resid   = None,
        cutoff         = cutoff,
    )


def build_tleap_resid_map(pdb_path: str,
                          chains: list = None) -> dict:
    """
    Build a mapping from (chain, original_resid) → tleap sequential resid.
    tleap renumbers residues sequentially across all chains.

    Parameters
    ----------
    pdb_path : str  — prepared PDB (protein only, no ligand)
    chains   : list — chains to include (default: all)

    Returns
    -------
    dict with keys:
        resid_map : {(chain, orig_resid): tleap_resid}
        warnings  : list
        errors    : list
    """
    result = {"resid_map": {}, "warnings": [], "errors": []}

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    seen    = []
    seen_set = set()

    for line in lines:
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        chain  = line[21]
        resseq = line[22:26].strip()
        icode  = line[26]
        key    = (chain, resseq, icode)

        if chains and chain not in chains:
            continue
        if key not in seen_set:
            seen.append((chain, resseq))
            seen_set.add(key)

    resid_map = {}
    for tleap_num, (chain, orig) in enumerate(seen, start=1):
        resid_map[(chain, orig)] = tleap_num

    result["resid_map"] = resid_map
    return result


def format_pacs_selection(binding_site_resids: list,
                          tleap_resid_map: dict = None) -> str:
    """
    Format binding site residues as PaCS-Q selection string.
    Uses tleap sequential resids if map is provided.

    Parameters
    ----------
    binding_site_resids : list of {chain, resid, resname, tleap_resid}
    tleap_resid_map     : {(chain, orig_resid): tleap_resid}

    Returns
    -------
    str — e.g. "resid 50 52 87 101"
    """
    resids = []
    for res in binding_site_resids:
        if res.get("tleap_resid") is not None:
            resids.append(int(res["tleap_resid"]))
        elif tleap_resid_map:
            key = (res["chain"], str(res["resid"]))
            if key in tleap_resid_map:
                resids.append(tleap_resid_map[key])
            else:
                resids.append(int(res["resid"]))
        else:
            resids.append(int(res["resid"]))

    resids_sorted = sorted(set(resids))
    return "resid " + " ".join(str(r) for r in resids_sorted)


def get_coords(pdb_path: str,
               record_types: tuple = ("ATOM", "HETATM"),
               skip_hydrogens: bool = True) -> np.ndarray:
    """
    Extract all heavy atom coordinates from a PDB file.

    Returns
    -------
    np.ndarray shape (N, 3)
    """
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in record_types:
                continue
            if skip_hydrogens:
                element  = line[76:78].strip().upper() if len(line) > 76 else ""
                atomname = line[12:16].strip()
                if element == "H" or atomname.startswith("H"):
                    continue
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except ValueError:
                continue
    return np.array(coords) if coords else np.zeros((0, 3))


def count_atoms(pdb_path: str,
                record_types: tuple = ("ATOM", "HETATM")) -> int:
    """Count heavy atoms in PDB file."""
    return len(get_coords(pdb_path, record_types=record_types))


def write_pdb(lines: list, output_path: str) -> None:
    """Write PDB lines to file, ensuring END record."""
    with open(output_path, "w") as f:
        f.writelines(lines)
        if lines and not lines[-1].strip().startswith("END"):
            f.write("END\n")


# ============================================================
# SECTION E — PROTONATION (PDB2PQR + CYX + HIS)
# ============================================================

def run_pdb2pqr(pdb_path: str, pH: float = 7.4,
                output_prefix: str = None) -> dict:
    """
    Run PDB2PQR with AMBER force field and propka for protonation.

    Parameters
    ----------
    pdb_path      : str   — input PDB
    pH            : float — target pH (default 7.4)
    output_prefix : str   — output prefix (default: stem_pH74)

    Returns
    -------
    dict with keys:
        pdb_path      : str  — protonated PDB path
        his_states    : list — [{resid, chain, state: HID/HIE/HIP}]
        warnings      : list
        errors        : list
    """
    result = {
        "pdb_path":   None,
        "his_states": [],
        "warnings":   [],
        "errors":     [],
    }

    if output_prefix is None:
        ph_tag        = str(pH).replace(".", "")
        output_prefix = f"{Path(pdb_path).stem}_pH{ph_tag}"

    out_pdb  = f"{output_prefix}.pdb"
    out_pqr  = f"{output_prefix}.pqr"
    log_file = f"{output_prefix}_pdb2pqr.log"

    cmd = [
        "pdb2pqr",
        "--ff=AMBER",
        f"--pH={pH}",
        "--titration-state-method=propka",
        "--with-ph=" + str(pH),
        "--pdb-output", out_pdb,
        pdb_path,
        out_pqr,
    ]

    try:
        with open(log_file, "w") as log:
            proc = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=120,
            )
        if proc.returncode != 0:
            result["errors"].append(
                f"PDB2PQR failed (return code {proc.returncode}). "
                f"Check log: {log_file}"
            )
            return result
    except FileNotFoundError:
        result["errors"].append(
            "pdb2pqr not found. Install with: pip install pdb2pqr"
        )
        return result
    except subprocess.TimeoutExpired:
        result["errors"].append("PDB2PQR timed out (>120s).")
        return result

    # Parse HIS protonation states from output PDB
    if os.path.exists(out_pdb):
        his_states = _parse_his_states(out_pdb)
        result["his_states"] = his_states
        result["pdb_path"]   = out_pdb
    else:
        result["errors"].append(
            f"PDB2PQR did not produce {out_pdb}."
        )

    return result


def _parse_his_states(pdb_path: str) -> list:
    """
    Parse HIS protonation states (HID/HIE/HIP) from PDB2PQR output.

    Returns
    -------
    list of {chain, resid, resname (HID/HIE/HIP)}
    """
    states = []
    seen   = set()

    with open(pdb_path, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec != "ATOM":
                continue
            resname = line[17:20].strip().upper()
            if resname not in ("HID", "HIE", "HIP"):
                continue
            chain  = line[21]
            resseq = line[22:26].strip()
            key    = (chain, resseq)
            if key not in seen:
                seen.add(key)
                states.append({
                    "chain":   chain,
                    "resid":   resseq,
                    "resname": resname,
                })

    return states


def override_his_protonation(pdb_path: str,
                              overrides: dict,
                              output_path: str = None) -> dict:
    """
    Override HIS protonation states in PDB.

    Parameters
    ----------
    pdb_path  : str  — input PDB (from PDB2PQR)
    overrides : dict — {(chain, resid): "HID" | "HIE" | "HIP"}
    output_path : str

    Returns
    -------
    dict with keys: pdb_path, warnings, errors
    """
    result = {"pdb_path": None, "warnings": [], "errors": []}

    if output_path is None:
        output_path = f"{Path(pdb_path).stem}_hisfix.pdb"

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        rec = line[:6].strip()
        if rec == "ATOM":
            resname = line[17:20].strip().upper()
            chain   = line[21]
            resseq  = line[22:26].strip()
            key     = (chain, resseq)
            if resname in ("HID", "HIE", "HIP") and key in overrides:
                new_res  = overrides[key].upper()
                new_line = line[:17] + f"{new_res:<3s}" + line[20:]
                new_lines.append(new_line)
                continue
        new_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(new_lines)

    result["pdb_path"] = output_path
    return result


def detect_disulfide_pairs(pdb_path: str) -> dict:
    """
    Detect disulfide bonds by measuring S-S distances between CYS residues.

    Thresholds:
        < 2.5 Å  → confirmed CYX pair
        2.5-3.0 Å → borderline — warn user, ask confirmation
        ≥ 3.0 Å  → not a disulfide

    Returns
    -------
    dict with keys:
        confirmed   : list of ((chainA, resA), (chainB, resB))
        borderline  : list of ((chainA, resA), (chainB, resB), distance)
        warnings    : list
        errors      : list
    """
    result = {
        "confirmed":  [],
        "borderline": [],
        "warnings":   [],
        "errors":     [],
    }

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Collect SG atoms of CYS/CYX residues
    sg_atoms = []
    for line in lines:
        rec = line[:6].strip()
        if rec != "ATOM":
            continue
        resname  = line[17:20].strip().upper()
        atomname = line[12:16].strip()
        if resname not in ("CYS", "CYX") or atomname != "SG":
            continue
        chain  = line[21]
        resseq = line[22:26].strip()
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            sg_atoms.append((chain, resseq, np.array([x, y, z])))
        except ValueError:
            continue

    # Pairwise distances
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            ch_a, res_a, coord_a = sg_atoms[i]
            ch_b, res_b, coord_b = sg_atoms[j]
            dist = float(np.linalg.norm(coord_b - coord_a))

            if dist < 2.5:
                result["confirmed"].append(
                    ((ch_a, res_a), (ch_b, res_b))
                )
            elif dist < 3.0:
                result["borderline"].append(
                    ((ch_a, res_a), (ch_b, res_b), round(dist, 2))
                )

    if result["borderline"]:
        bl = [(f"{a[0]}:{a[1]}-{b[0]}:{b[1]} ({d:.2f}Å)")
              for a, b, d in result["borderline"]]
        result["warnings"].append(
            f"Borderline S-S distances found: {bl}. "
            "Please confirm if these are disulfide bonds."
        )

    return result


# ============================================================
# SECTION F — COFACTOR HANDLING
# ============================================================

# Pre-built AMBER parameter URLs for common cofactors
COFACTOR_PARAM_URLS = {
    "NAD": {
        "prep":   "http://amber.manchester.ac.uk/amber/NAD/NAD.prep",
        "frcmod": "http://amber.manchester.ac.uk/amber/NAD/NAD.frcmod",
    },
    "FAD": {
        "prep":   "http://amber.manchester.ac.uk/amber/FAD/FAD.prep",
        "frcmod": "http://amber.manchester.ac.uk/amber/FAD/FAD.frcmod",
    },
    # Add more as needed
}


def detect_cofactors(pdb_path: str,
                     hetatm_classification: dict = None) -> dict:
    """
    Detect known cofactors and classify them for parameter handling.

    Returns
    -------
    dict with keys:
        cofactors       : list of {resname, chain, resid, param_status}
            param_status: "known_prebuilt" | "unknown_need_params" | "gaff2_fallback"
        warnings        : list
        errors          : list
    """
    result = {"cofactors": [], "warnings": [], "errors": []}

    if hetatm_classification is None:
        hetatm_classification = classify_hetatm(pdb_path)

    for entry in hetatm_classification.get("cofactors", []):
        resname = entry["resname"]
        if resname in COFACTOR_PARAM_URLS:
            status = "known_prebuilt"
        else:
            status = "unknown_need_params"
            result["warnings"].append(
                f"Unknown cofactor {resname} — parameters not in library. "
                "User can upload .prep+.frcmod or GAFF2 will be used as fallback."
            )
        result["cofactors"].append({**entry, "param_status": status})

    return result


def download_amber_params(resname: str,
                          output_dir: str = ".") -> dict:
    """
    Download prebuilt AMBER parameters for known cofactors.

    Returns
    -------
    dict with keys: prep_path, frcmod_path, warnings, errors
    """
    result = {
        "prep_path":   None,
        "frcmod_path": None,
        "warnings":    [],
        "errors":      [],
    }

    resname = resname.upper()
    if resname not in COFACTOR_PARAM_URLS:
        result["errors"].append(
            f"No prebuilt parameters for {resname}. "
            "Use GAFF2 fallback or upload manually."
        )
        return result

    urls    = COFACTOR_PARAM_URLS[resname]
    prep_p  = os.path.join(output_dir, f"{resname}.prep")
    frcm_p  = os.path.join(output_dir, f"{resname}.frcmod")

    for url, dest in [(urls["prep"], prep_p), (urls["frcmod"], frcm_p)]:
        try:
            urllib.request.urlretrieve(url, dest)
        except Exception as e:
            result["errors"].append(f"Download failed for {url}: {e}")
            return result

    result["prep_path"]   = prep_p
    result["frcmod_path"] = frcm_p
    return result


def check_atom_names(pdb_path: str,
                     prep_path: str,
                     resname: str) -> dict:
    """
    Check that atom names in PDB match those in the .prep file.
    Mismatches cause tleap to fail silently or crash.

    Returns
    -------
    dict with keys:
        matched   : bool
        missing   : list — atoms in prep but not in PDB
        extra     : list — atoms in PDB but not in prep
        warnings  : list
        errors    : list
    """
    result = {
        "matched":  False,
        "missing":  [],
        "extra":    [],
        "warnings": [],
        "errors":   [],
    }

    # Collect atom names from PDB
    pdb_atoms = set()
    with open(pdb_path, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip().upper() != resname.upper():
                continue
            pdb_atoms.add(line[12:16].strip())

    # Collect atom names from prep
    prep_atoms = set()
    with open(prep_path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                prep_atoms.add(parts[1])

    missing = list(prep_atoms - pdb_atoms)
    extra   = list(pdb_atoms - prep_atoms)

    if missing or extra:
        result["warnings"].append(
            f"Atom name mismatch for {resname}: "
            f"missing={missing}, extra={extra}"
        )
    else:
        result["matched"] = True

    result["missing"] = missing
    result["extra"]   = extra
    return result


def auto_fix_atom_names(pdb_path: str,
                        resname: str,
                        name_map: dict,
                        output_path: str = None) -> dict:
    """
    Rename atoms in PDB according to name_map.

    Parameters
    ----------
    name_map : {old_name: new_name}
    """
    result = {"pdb_path": None, "warnings": [], "errors": []}

    if output_path is None:
        output_path = f"{Path(pdb_path).stem}_fixnames.pdb"

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            if line[17:20].strip().upper() == resname.upper():
                old_name = line[12:16].strip()
                if old_name in name_map:
                    new_name = f"{name_map[old_name]:<4s}"
                    line     = line[:12] + new_name + line[16:]
        new_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(new_lines)

    result["pdb_path"] = output_path
    return result


# ============================================================
# SECTION G — LIGAND PREPARATION (4 MODES)
# ============================================================

def force_resname_mol(file_path: str, file_type: str = "pdb") -> dict:
    """
    Rename all ligand residue names to MOL in-place.
    Applies to PDB, prep, mol2, and lib files.

    Returns
    -------
    dict with keys: file_path, warnings, errors
    """
    result = {"file_path": file_path, "warnings": [], "errors": []}

    with open(file_path, "r") as f:
        content = f.read()

    file_type = file_type.lower()

    if file_type == "pdb":
        lines     = content.splitlines(keepends=True)
        new_lines = []
        for line in lines:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                line = line[:17] + "MOL" + line[20:]
            new_lines.append(line)
        new_content = "".join(new_lines)

    elif file_type == "prep":
        # Replace residue name in PREP header section
        lines     = content.splitlines(keepends=True)
        new_lines = []
        header_done = False
        for line in lines:
            if not header_done and len(line) > 3:
                # Prep files have resname in 3rd column of atom table header
                # and in the 2nd line after blank
                stripped = line.strip()
                if stripped and len(stripped.split()) >= 1:
                    # Heuristic: short 3-letter word alone on a line = resname
                    parts = stripped.split()
                    if (len(parts) == 1 and len(parts[0]) <= 4
                            and parts[0].isalpha()):
                        new_lines.append("MOL\n")
                        header_done = True
                        continue
            new_lines.append(line)
        new_content = "".join(new_lines)

    elif file_type == "mol2":
        # Replace @<TRIPOS>MOLECULE section name
        lines     = content.splitlines(keepends=True)
        new_lines = []
        in_mol_section = False
        mol_name_written = False
        for line in lines:
            if "@<TRIPOS>MOLECULE" in line:
                in_mol_section  = True
                mol_name_written = False
                new_lines.append(line)
                continue
            if in_mol_section and not mol_name_written:
                new_lines.append("MOL\n")
                mol_name_written = True
                in_mol_section   = False
                continue
            # Also fix ATOM section residue names (col 8 in mol2)
            if "@<TRIPOS>ATOM" in line:
                new_lines.append(line)
                continue
            new_lines.append(line)
        new_content = "".join(new_lines)

    else:
        result["warnings"].append(
            f"force_resname_mol: unsupported file type '{file_type}'. Skipping."
        )
        return result

    with open(file_path, "w") as f:
        f.write(new_content)

    return result


def validate_prep_file(prep_path: str) -> dict:
    """
    Validate that a .prep file has required sections:
    LOOP, DONE, and atom table.

    Returns
    -------
    dict with keys: valid, warnings, errors
    """
    result = {"valid": False, "warnings": [], "errors": []}

    with open(prep_path, "r") as f:
        content = f.read()

    if "LOOP" not in content:
        result["errors"].append(
            f"{prep_path}: Missing LOOP section. "
            "File may be corrupt or incomplete."
        )
    if "DONE" not in content:
        result["errors"].append(
            f"{prep_path}: Missing DONE section."
        )
    if not result["errors"]:
        result["valid"] = True

    return result


def prep_to_pdb(prep_path: str, output_pdb: str = None) -> dict:
    """
    Convert .prep file to PDB using antechamber.

    Returns
    -------
    dict with keys: pdb_path, warnings, errors
    """
    result = {"pdb_path": None, "warnings": [], "errors": []}

    if output_pdb is None:
        output_pdb = f"{Path(prep_path).stem}_from_prep.pdb"

    cmd = [
        "antechamber",
        "-i",  prep_path, "-fi", "prepi",
        "-o",  output_pdb, "-fo", "pdb",
        "-dr", "no",
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if proc.returncode != 0 or not os.path.exists(output_pdb):
            result["errors"].append(
                f"antechamber prep→pdb failed: {proc.stderr[:300]}"
            )
            return result
    except FileNotFoundError:
        result["errors"].append("antechamber not found. Ensure AmberTools is installed.")
        return result

    result["pdb_path"] = output_pdb
    return result


def mol2_has_charges(mol2_path: str) -> bool:
    """
    Check if a .mol2 file already has partial charges
    (non-zero values in charge column of ATOM section).
    """
    in_atom = False
    with open(mol2_path, "r") as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if "@<TRIPOS>" in line and in_atom:
                break
            if in_atom:
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charge = float(parts[8])
                        if abs(charge) > 1e-6:
                            return True
                    except ValueError:
                        pass
    return False


def lib_needs_frcmod(lib_path: str) -> bool:
    """
    Check if a .lib file references non-standard GAFF2 atom types
    that would require a separate .frcmod.
    Returns True if frcmod is likely needed.
    """
    with open(lib_path, "r") as f:
        content = f.read()
    # Heuristic: if atom types beyond standard GAFF2 set appear
    # (we check for types not in standard set)
    standard_gaff2_types = {
        "c", "c1", "c2", "c3", "ca", "cp", "cq", "cc", "cd",
        "ce", "cf", "cg", "ch", "cx", "cy", "cu", "cv", "cz",
        "h1", "h2", "h3", "h4", "h5", "ha", "hc", "hn", "ho",
        "hp", "hs", "hw", "hx", "f", "cl", "br", "i",
        "n", "n1", "n2", "n3", "n4", "na", "nb", "nc", "nd",
        "ne", "nf", "nh", "no", "ns", "nt", "nu", "nv", "nx",
        "ny", "nz", "n+", "n8", "o", "oh", "os", "ow", "p2",
        "p3", "p4", "p5", "pb", "pc", "pd", "pe", "pf", "px",
        "py", "s", "s2", "s4", "s6", "sh", "ss", "sx", "sy",
    }
    # Extract atom type column from !entry section
    for line in content.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].startswith('"'):
            atype = parts[-1].strip('"').lower()
            if atype and atype not in standard_gaff2_types:
                return True
    return False


def validate_charge(total_charge: float,
                    threshold: float = 4.0) -> dict:
    """
    Warn if absolute charge exceeds threshold.
    Most drug-like molecules have |charge| ≤ 2.
    """
    result = {"warnings": [], "errors": []}
    if abs(total_charge) > threshold:
        result["warnings"].append(
            f"Ligand charge = {total_charge:.1f}. "
            f"|charge| > {threshold} is unusual for drug-like molecules. "
            "Please verify the protonation state."
        )
    return result


def check_ph_consistency(protein_pH: float,
                         ligand_pH: float,
                         tolerance: float = 0.5) -> dict:
    """
    Check that protein and ligand were prepared at the same pH.
    """
    result = {"consistent": True, "warnings": [], "errors": []}
    if abs(protein_pH - ligand_pH) > tolerance:
        result["consistent"] = False
        result["warnings"].append(
            f"pH mismatch: protein prepared at pH {protein_pH}, "
            f"ligand at pH {ligand_pH}. "
            "Use the same pH for consistent protonation states."
        )
    return result


def optimize_ligand_geometry(pdb_path: str,
                              output_path: str = None) -> dict:
    """
    Optimize ligand geometry using MMFF94 force field via RDKit.
    Used after prep→pdb or mol2→pdb conversion.

    Returns
    -------
    dict with keys: pdb_path, warnings, errors
    """
    result = {"pdb_path": None, "warnings": [], "errors": []}

    if output_path is None:
        output_path = f"{Path(pdb_path).stem}_opt.pdb"

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol is None:
            result["errors"].append(
                f"RDKit could not read {pdb_path}."
            )
            return result

        mol = Chem.AddHs(mol, addCoords=True)
        ff  = AllChem.MMFFGetMoleculeForceField(
            mol, AllChem.MMFFGetMoleculeProperties(mol)
        )
        if ff is None:
            result["warnings"].append(
                "MMFF94 not available for this molecule. "
                "Skipping geometry optimization."
            )
            shutil.copy(pdb_path, output_path)
            result["pdb_path"] = output_path
            return result

        ff.Minimize()
        Chem.MolToPDBFile(mol, output_path)

    except ImportError:
        result["warnings"].append(
            "RDKit not available. Skipping geometry optimization."
        )
        shutil.copy(pdb_path, output_path)

    result["pdb_path"] = output_path
    return result


def smiles_to_amber(smiles: str,
                    pH:    float = 7.4,
                    charge: int  = None,
                    name:   str  = "guest") -> dict:
    """
    MODE 1: SMILES → GAFF2 parameters (AM1-BCC charges).

    Steps:
        1. Protonate with dimorphite-dl at given pH
        2. Generate 3D structure with RDKit ETKDG + UFF
        3. antechamber AM1-BCC → GAFF2 → .prep + .frcmod
        4. force_resname_mol() → rename to MOL

    Returns
    -------
    dict with keys:
        prep_path    : str
        frcmod_path  : str
        pdb_path     : str  — 3D geometry for BSV
        charge       : int
        warnings     : list
        errors       : list
    """
    result = {
        "prep_path":   None,
        "frcmod_path": None,
        "pdb_path":    None,
        "charge":      None,
        "warnings":    [],
        "errors":      [],
    }

    smiles = smiles.strip()

    # Step 1: Protonate at pH
    try:
        from dimorphite_dl import protonate_smiles
        protonated = protonate_smiles(smiles, min_ph=pH - 0.5,
                                      max_ph=pH + 0.5)
        if protonated:
            smiles = protonated[0]
            if len(protonated) > 1:
                result["warnings"].append(
                    f"Multiple protonation states at pH {pH}: "
                    f"{protonated}. Using first state."
                )
    except ImportError:
        result["warnings"].append(
            "dimorphite-dl not available. Using input SMILES as-is."
        )

    # Step 2: 3D structure
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["errors"].append("Invalid SMILES — RDKit could not parse.")
            return result

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        if charge is None:
            charge = Chem.GetFormalCharge(mol)

        mol_file = f"{name}.mol"
        pdb_file = f"{name}.pdb"
        Chem.MolToMolFile(mol, mol_file)
        Chem.MolToPDBFile(mol, pdb_file)

    except ImportError:
        result["errors"].append("RDKit not available.")
        return result

    # Step 3: antechamber AM1-BCC
    prep_file   = f"{name}.prep"
    frcmod_file = f"{name}.frcmod"
    log_file    = f"{name}_antechamber.log"

    cmd = [
        "antechamber",
        "-i",  mol_file, "-fi", "mdl",
        "-o",  prep_file, "-fo", "prepi",
        "-rn", "MOL",
        "-at", "gaff2",
        "-c",  "bcc",
        "-nc", str(charge),
        "-s",  "2",
    ]

    try:
        with open(log_file, "w") as log:
            proc = subprocess.run(
                cmd, stdout=log, stderr=subprocess.STDOUT,
                timeout=1200,
            )
        if proc.returncode != 0:
            result["errors"].append(
                f"antechamber AM1-BCC failed. Check log: {log_file}"
            )
            return result
    except FileNotFoundError:
        result["errors"].append("antechamber not found.")
        return result
    except subprocess.TimeoutExpired:
        result["errors"].append("antechamber timed out (>1200s).")
        return result

    # parmchk2
    cmd2 = ["parmchk2", "-i", prep_file, "-f", "prepi", "-o", frcmod_file]
    try:
        subprocess.run(cmd2, check=True, timeout=60)
    except Exception as e:
        result["errors"].append(f"parmchk2 failed: {e}")
        return result

    result["prep_path"]   = prep_file
    result["frcmod_path"] = frcmod_file
    result["pdb_path"]    = pdb_file
    result["charge"]      = charge
    return result


def prep_frcmod_to_amber(prep_path: str,
                         frcmod_path: str,
                         pH: float = 7.4) -> dict:
    """
    MODE 2: .prep + .frcmod upload → validate + generate PDB.

    Returns
    -------
    dict with keys: prep_path, frcmod_path, pdb_path, warnings, errors
    """
    result = {
        "prep_path":   prep_path,
        "frcmod_path": frcmod_path,
        "pdb_path":    None,
        "warnings":    [],
        "errors":      [],
    }

    # Validate
    val = validate_prep_file(prep_path)
    result["errors"]   += val["errors"]
    result["warnings"] += val["warnings"]
    if val["errors"]:
        return result

    # Convert prep → PDB
    p2p = prep_to_pdb(prep_path)
    result["errors"]   += p2p["errors"]
    result["warnings"] += p2p["warnings"]
    if p2p["errors"]:
        return result

    # Optimize geometry
    opt = optimize_ligand_geometry(p2p["pdb_path"])
    result["warnings"] += opt["warnings"]
    result["errors"]   += opt["errors"]

    pdb_path = opt["pdb_path"] or p2p["pdb_path"]

    # Rename to MOL
    force_resname_mol(prep_path,  "prep")
    force_resname_mol(pdb_path,   "pdb")

    result["pdb_path"] = pdb_path
    return result


def mol2_to_amber(mol2_path: str,
                  pH: float = 7.4,
                  charge: int = None) -> dict:
    """
    MODE 3: .mol2 → GAFF2 parameters.

    Returns
    -------
    dict with keys: prep_path, frcmod_path, pdb_path, charge, warnings, errors
    """
    result = {
        "prep_path":   None,
        "frcmod_path": None,
        "pdb_path":    None,
        "charge":      charge,
        "warnings":    [],
        "errors":      [],
    }

    has_charges = mol2_has_charges(mol2_path)
    stem        = Path(mol2_path).stem
    prep_file   = f"{stem}.prep"
    frcmod_file = f"{stem}.frcmod"
    pdb_file    = f"{stem}.pdb"

    charge_method = "wc" if has_charges else "bcc"
    if not has_charges:
        result["warnings"].append(
            "mol2 file has no partial charges — running AM1-BCC."
        )

    nc_str = str(charge) if charge is not None else "0"

    cmd = [
        "antechamber",
        "-i",  mol2_path, "-fi", "mol2",
        "-o",  prep_file, "-fo", "prepi",
        "-rn", "MOL",
        "-at", "gaff2",
        "-c",  charge_method,
        "-nc", nc_str,
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True,
                              text=True, timeout=1200)
        if proc.returncode != 0:
            result["errors"].append(
                f"antechamber mol2 failed: {proc.stderr[:300]}"
            )
            return result
    except Exception as e:
        result["errors"].append(f"antechamber error: {e}")
        return result

    # parmchk2
    cmd2 = ["parmchk2", "-i", prep_file, "-f", "prepi", "-o", frcmod_file]
    try:
        subprocess.run(cmd2, check=True, timeout=60)
    except Exception as e:
        result["errors"].append(f"parmchk2 failed: {e}")
        return result

    # mol2 → pdb via obabel
    cmd3 = ["obabel", mol2_path, "-O", pdb_file]
    try:
        subprocess.run(cmd3, check=True, timeout=30)
    except Exception as e:
        result["errors"].append(f"obabel mol2→pdb failed: {e}")
        return result

    force_resname_mol(pdb_file, "pdb")
    force_resname_mol(prep_file, "prep")

    result["prep_path"]   = prep_file
    result["frcmod_path"] = frcmod_file
    result["pdb_path"]    = pdb_file
    return result


def lib_to_amber(lib_path: str,
                 frcmod_path: str = None) -> dict:
    """
    MODE 4: .lib (AMBER library) → generate PDB for BSV.

    Returns
    -------
    dict with keys: lib_path, frcmod_path, pdb_path, warnings, errors
    """
    result = {
        "lib_path":    lib_path,
        "frcmod_path": frcmod_path,
        "pdb_path":    None,
        "warnings":    [],
        "errors":      [],
    }

    stem     = Path(lib_path).stem
    pdb_file = f"{stem}.pdb"

    needs_frc = lib_needs_frcmod(lib_path)
    if needs_frc and frcmod_path is None:
        result["warnings"].append(
            "This .lib file may contain non-standard atom types. "
            "A .frcmod file is recommended. "
            "Proceeding without it — tleap may fail."
        )

    # lib → pdb via obabel
    cmd = ["obabel", lib_path, "-O", pdb_file]
    try:
        subprocess.run(cmd, check=True, timeout=30)
    except Exception as e:
        result["errors"].append(f"obabel lib→pdb failed: {e}")
        return result

    force_resname_mol(pdb_file, "pdb")

    result["pdb_path"] = pdb_file
    return result


def generate_ligand_3d_pdb(source: str,
                            file_type: str,
                            pH: float = 7.4) -> dict:
    """
    Unified entry point for ligand 3D PDB generation.
    Dispatches to correct mode based on file_type.

    Parameters
    ----------
    source    : str  — SMILES string or file path
    file_type : str  — "smiles" | "prep" | "mol2" | "lib"
    pH        : float

    Returns
    -------
    dict with pdb_path + mode-specific keys
    """
    file_type = file_type.lower()

    if file_type == "smiles":
        return smiles_to_amber(source, pH=pH)
    elif file_type == "prep":
        # Expect source to be "prep_path,frcmod_path"
        parts = source.split(",")
        if len(parts) != 2:
            return {"errors": ["For prep mode, provide 'prep_path,frcmod_path'"],
                    "warnings": []}
        return prep_frcmod_to_amber(parts[0].strip(), parts[1].strip(), pH)
    elif file_type == "mol2":
        return mol2_to_amber(source, pH=pH)
    elif file_type == "lib":
        return lib_to_amber(source)
    else:
        return {
            "errors":   [f"Unknown file_type: {file_type}"],
            "warnings": [],
        }


# ============================================================
# SECTION H — BSV PLACEMENT GEOMETRY
# ============================================================

def _normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector."""
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else v


def _subdivide_icosphere(vertices: np.ndarray,
                          faces: np.ndarray,
                          subdivisions: int) -> tuple:
    """
    Subdivide icosphere faces by splitting each triangle edge
    and projecting new vertices onto the unit sphere.

    Returns
    -------
    vertices : np.ndarray (N, 3)
    faces    : np.ndarray (M, 3)
    """
    for _ in range(subdivisions):
        new_faces    = []
        midpoint_cache = {}


        def midpoint(i, j):
            nonlocal vertices                      # ← FIRST line
            key = (min(i, j), max(i, j))
            if key in midpoint_cache:
                return midpoint_cache[key]
            mid      = _normalize((vertices[i] + vertices[j]) / 2)
            vertices = np.vstack([vertices, mid])
            idx      = len(vertices) - 1
            midpoint_cache[key] = idx
            return idx
      

        for tri in faces:
            a, b, c = tri
            ab = midpoint(a, b)
            bc = midpoint(b, c)
            ca = midpoint(c, a)
            new_faces.extend([
                [a, ab, ca],
                [b, bc, ab],
                [c, ca, bc],
                [ab, bc, ca],
            ])
        faces = np.array(new_faces)

    return vertices, faces


def generate_icosphere(subdivisions: int = 3) -> np.ndarray:
    """
    Generate unit icosphere vertices (triangular lattice on unit sphere).
    Each vertex defines a vector direction for BSV.

    Parameters
    ----------
    subdivisions : int — number of face subdivision steps
        0 → 12 vertices (icosahedron)
        1 → 42 vertices
        2 → 162 vertices
        3 → 642 vertices  ← default (≈ 492 usable before hemisphere filter)

    Returns
    -------
    np.ndarray shape (N, 3) — unit vectors on sphere surface
    """
    # Golden ratio
    phi = (1.0 + math.sqrt(5.0)) / 2.0

    # 12 icosahedron vertices
    verts = np.array([
        [-1,  phi, 0], [ 1,  phi, 0], [-1, -phi, 0], [ 1, -phi, 0],
        [ 0, -1,  phi], [ 0,  1,  phi], [ 0, -1, -phi], [ 0,  1, -phi],
        [ phi, 0, -1], [ phi, 0,  1], [-phi, 0, -1], [-phi, 0,  1],
    ], dtype=float)

    # Normalize to unit sphere
    verts = np.array([_normalize(v) for v in verts])

    # 20 icosahedron faces
    faces = np.array([
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ], dtype=int)

    verts, faces = _subdivide_icosphere(verts, faces, subdivisions)

    # Return unique vertices (deduplicate)
    verts_rounded = np.round(verts, 8)
    _, unique_idx = np.unique(verts_rounded, axis=0, return_index=True)
    return verts[unique_idx]


def _get_vdw_radius(element: str) -> float:
    """Return van der Waals radius for element (Å). Default 1.70 if unknown."""
    return VDW_RADII.get(element.upper().strip(), 1.70)


def _bsv_intersect(origin: np.ndarray,
                   direction: np.ndarray,
                   atom_center: np.ndarray,
                   atom_vdw: float) -> Optional[float]:
    """
    Compute BSV vector length c (intersection with atom vdW surface).

    Geometry (Kuvek et al. Fig. 3):
        a = |origin → atom_center|
        b = a · cos(α)  where α = angle between v_black and v_blue
        y = √(a² - b²)  perpendicular distance
        x = √(r² - y²)  if y ≤ r (intersection exists)
        c = b - x        final vector length

    Returns
    -------
    float c if intersection exists, else None
    """
    v_black = atom_center - origin
    a       = float(np.linalg.norm(v_black))

    if a < 1e-10:
        return None

    cos_alpha = float(np.dot(v_black, direction) / a)
    b         = a * cos_alpha

    # b must be > 0 (atom is in front of origin along direction)
    if b <= 0:
        return None

    y_sq = a * a - b * b
    r    = atom_vdw

    if y_sq > r * r:
        return None  # vector misses atom

    y = math.sqrt(max(y_sq, 0.0))
    x = math.sqrt(max(r * r - y * y, 0.0))
    c = b - x

    return c if c > 0 else None


def compute_bsv_vectors(origin:      np.ndarray,
                         protein_pdb: str,
                         sphere_r:    float = 20.0,
                         subdivisions: int  = 3,
                         hemisphere:  bool  = True) -> dict:
    """
    Compute Binding Site Vectors anchored at `origin`.
    Implements Section 2.1 of Kuvek et al. JCIM 2026.

    Parameters
    ----------
    origin       : np.ndarray (3,) — anchor point (e.g. binding site center)
    protein_pdb  : str             — protein PDB (no ligand)
    sphere_r     : float           — initial sphere radius in Å
    subdivisions : int             — icosphere subdivisions
    hemisphere   : bool            — restrict to +z hemisphere (above heme analog)

    Returns
    -------
    dict with keys:
        vectors     : np.ndarray (N, 3) — direction unit vectors
        lengths     : np.ndarray (N,)   — final vector lengths (Å)
        charges     : np.ndarray (N,)   — partial charges at intersection
        n_vectors   : int
        warnings    : list
        errors      : list
    """
    result = {
        "vectors":   None,
        "lengths":   None,
        "charges":   None,
        "n_vectors": 0,
        "warnings":  [],
        "errors":    [],
    }

    # Generate icosphere directions
    directions = generate_icosphere(subdivisions)

    # Optionally restrict to hemisphere (z > 0 relative to origin)
    if hemisphere:
        directions = directions[directions[:, 2] >= 0]

    n_vec = len(directions)

    # Collect protein heavy atoms: (coord, vdw_radius, partial_charge)
    atoms = []
    with open(protein_pdb, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            element  = line[76:78].strip().upper() if len(line) > 76 else ""
            atomname = line[12:16].strip()
            if element == "H" or atomname.startswith("H"):
                continue
            if not element:
                # Guess from atom name
                element = re.sub(r'[^A-Z]', '', atomname[:2])

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            # Partial charge: B-factor column as placeholder
            # (real charges assigned by GROMOS/AMBER FF in full impl)
            try:
                charge = float(line[60:66].strip() or "0.0")
            except (ValueError, IndexError):
                charge = 0.0

            vdw = _get_vdw_radius(element)
            atoms.append((np.array([x, y, z]), vdw, charge))

    if not atoms:
        result["errors"].append("No protein atoms found in PDB.")
        return result

    # For each vector direction, find closest intersection
    lengths  = np.full(n_vec, sphere_r)   # default: sphere radius
    charges  = np.zeros(n_vec)

    for k, direction in enumerate(directions):
        best_c      = sphere_r
        best_charge = 0.0

        for (atom_coord, atom_vdw, atom_charge) in atoms:
            c = _bsv_intersect(origin, direction, atom_coord, atom_vdw)
            if c is not None and c < best_c:
                best_c      = c
                best_charge = atom_charge

        lengths[k]  = best_c
        charges[k]  = best_charge

    result["vectors"]   = directions
    result["lengths"]   = lengths
    result["charges"]   = charges
    result["n_vectors"] = n_vec
    return result


def find_exit_direction(bsv_result:  dict,
                        ligand_pdb:  str,
                        origin:      np.ndarray) -> dict:
    """
    Find the BSV direction that points most toward the co-crystal ligand centroid.
    This is the natural exit direction of the pocket.

    Parameters
    ----------
    bsv_result  : dict — output of compute_bsv_vectors()
    ligand_pdb  : str  — co-crystal ligand PDB (heavy atoms)
    origin      : np.ndarray (3,)

    Returns
    -------
    dict with keys:
        exit_direction : np.ndarray (3,) — unit vector
        exit_length    : float           — BSV length in exit direction (Å)
        exit_index     : int             — vector index
        ligand_centroid: np.ndarray (3,)
        warnings       : list
        errors         : list
    """
    result = {
        "exit_direction":  None,
        "exit_length":     None,
        "exit_index":      None,
        "ligand_centroid": None,
        "warnings":        [],
        "errors":          [],
    }

    # Ligand centroid
    lig_coords = get_coords(ligand_pdb, record_types=("ATOM", "HETATM"))
    if len(lig_coords) == 0:
        result["errors"].append("No atoms found in ligand PDB.")
        return result

    centroid        = np.mean(lig_coords, axis=0)
    to_ligand       = _normalize(centroid - origin)

    vectors = bsv_result["vectors"]
    lengths = bsv_result["lengths"]

    # Find BSV closest in direction to ligand centroid
    dots     = vectors @ to_ligand
    best_idx = int(np.argmax(dots))

    result["exit_direction"]  = vectors[best_idx]
    result["exit_length"]     = float(lengths[best_idx])
    result["exit_index"]      = best_idx
    result["ligand_centroid"] = centroid
    return result


def get_ligand_radius(ligand_pdb: str) -> float:
    """
    Estimate ligand radius as max distance from centroid to any heavy atom.

    Returns
    -------
    float — radius in Å
    """
    coords = get_coords(ligand_pdb, record_types=("ATOM", "HETATM"))
    if len(coords) == 0:
        return 5.0  # safe default
    centroid = np.mean(coords, axis=0)
    dists    = np.linalg.norm(coords - centroid, axis=1)
    return float(np.max(dists))


def calculate_placement_distance(protein_pdb:  str,
                                  ligand_pdb:   str,
                                  exit_dir:     np.ndarray,
                                  origin:       np.ndarray,
                                  mode:         str = "protein_surface",
                                  buffer:       float = None,
                                  ref_resids:   list  = None,
                                  exit_length:  float = None) -> float:
    """
    Calculate final ligand placement distance from origin.

    Modes
    -----
    "protein_surface" (default):
        distance = exit_length + ligand_radius + buffer
    "specific_residues":
        distance = |origin → ref_resids_centroid| + buffer

    Parameters
    ----------
    buffer : float — placement buffer in Å (default: SILENT_DEFAULTS["default_buffer"])

    Returns
    -------
    float — distance from origin in Å (minimum: SILENT_DEFAULTS["min_placement_dist"])
    """
    if buffer is None:
        buffer = SILENT_DEFAULTS["default_buffer"]

    min_dist = SILENT_DEFAULTS["min_placement_dist"]

    if mode == "protein_surface":
        if exit_length is None:
            # Fallback: compute from protein surface
            bsv = compute_bsv_vectors(origin, protein_pdb,
                                      hemisphere=False)
            # Find max extent in exit direction
            dots = bsv["vectors"] @ exit_dir
            idx  = int(np.argmax(dots))
            exit_length = float(bsv["lengths"][idx])

        lig_radius = get_ligand_radius(ligand_pdb)
        distance   = exit_length + lig_radius + buffer

    elif mode == "specific_residues" and ref_resids:
        # Centroid of reference residue Cα atoms
        coords = get_coords(protein_pdb, record_types=("ATOM",))
        # Use all protein Cα as proxy (full implementation filters by resid)
        if len(coords) > 0:
            ref_centroid = np.mean(coords, axis=0)
            distance     = float(np.linalg.norm(ref_centroid - origin)) + buffer
        else:
            distance = buffer

    else:
        distance = buffer

    return max(distance, min_dist)


def place_ligand(ligand_pdb:   str,
                 origin:       np.ndarray,
                 exit_dir:     np.ndarray,
                 distance:     float,
                 output_path:  str = "ligand_placed.pdb") -> dict:
    """
    Translate ligand centroid to: origin + exit_dir × distance.

    Parameters
    ----------
    ligand_pdb  : str          — input ligand PDB
    origin      : np.ndarray   — BSV anchor point
    exit_dir    : np.ndarray   — unit exit direction
    distance    : float        — placement distance from origin (Å)
    output_path : str

    Returns
    -------
    dict with keys: pdb_path, placement_point, warnings, errors
    """
    result = {
        "pdb_path":        None,
        "placement_point": None,
        "warnings":        [],
        "errors":          [],
    }

    coords = get_coords(ligand_pdb, record_types=("ATOM", "HETATM"),
                        skip_hydrogens=False)
    if len(coords) == 0:
        result["errors"].append("No atoms in ligand PDB.")
        return result

    centroid        = np.mean(get_coords(ligand_pdb, record_types=("ATOM","HETATM")),
                              axis=0)
    target_position = origin + exit_dir * distance
    translation     = target_position - centroid

    # Apply translation to all atoms
    with open(ligand_pdb, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            try:
                x = float(line[30:38]) + translation[0]
                y = float(line[38:46]) + translation[1]
                z = float(line[46:54]) + translation[2]
                new_line = (line[:30]
                            + f"{x:8.3f}{y:8.3f}{z:8.3f}"
                            + line[54:])
                new_lines.append(new_line)
            except ValueError:
                new_lines.append(line)
        else:
            new_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(new_lines)
        if not new_lines[-1].strip().startswith("END"):
            f.write("END\n")

    result["pdb_path"]        = output_path
    result["placement_point"] = target_position
    return result


def check_clashes(protein_pdb: str,
                  ligand_pdb:  str,
                  threshold:   float = 1.5) -> dict:
    """
    Check for steric clashes between protein and placed ligand.
    A clash is defined as any protein-ligand heavy atom pair
    closer than `threshold` Å.

    Returns
    -------
    dict with keys:
        has_clashes  : bool
        n_clashes    : int
        clash_pairs  : list of (protein_atom_idx, ligand_atom_idx, distance)
        warnings     : list
        errors       : list
    """
    result = {
        "has_clashes": False,
        "n_clashes":   0,
        "clash_pairs": [],
        "warnings":    [],
        "errors":      [],
    }

    prot_coords = get_coords(protein_pdb, record_types=("ATOM",))
    lig_coords  = get_coords(ligand_pdb,  record_types=("ATOM", "HETATM"))

    if len(prot_coords) == 0 or len(lig_coords) == 0:
        return result

    # Vectorized pairwise distances
    for j, lc in enumerate(lig_coords):
        dists = np.linalg.norm(prot_coords - lc, axis=1)
        clash_idx = np.where(dists < threshold)[0]
        for i in clash_idx:
            result["clash_pairs"].append((int(i), j, float(dists[i])))

    if result["clash_pairs"]:
        result["has_clashes"] = True
        result["n_clashes"]   = len(result["clash_pairs"])
        result["warnings"].append(
            f"{result['n_clashes']} steric clashes detected "
            f"(threshold {threshold} Å). "
            "Buffer will be increased automatically."
        )

    return result


# ============================================================
# SECTION I — COMPLEX CENTERING + TLEAP
# ============================================================

def center_complex_at_origin(protein_pdb: str,
                              ligand_pdb:  str,
                              output_path: str = "complex_centered.pdb") -> dict:
    """
    Translate protein + ligand as a unit so their combined centroid
    is at the coordinate origin (0, 0, 0).
    Writes a single combined PDB file.

    Returns
    -------
    dict with keys: pdb_path, centroid, warnings, errors
    """
    result = {
        "pdb_path":  None,
        "centroid":  None,
        "warnings":  [],
        "errors":    [],
    }

    prot_coords = get_coords(protein_pdb, record_types=("ATOM",))
    lig_coords  = get_coords(ligand_pdb,  record_types=("ATOM", "HETATM"))

    if len(prot_coords) == 0:
        result["errors"].append("No protein atoms found.")
        return result

    all_coords = (np.vstack([prot_coords, lig_coords])
                  if len(lig_coords) > 0 else prot_coords)
    centroid   = np.mean(all_coords, axis=0)
    translation = -centroid

    def _translate_pdb(in_path, translation, record_types, out_lines):
        with open(in_path, "r") as f:
            for line in f:
                rec = line[:6].strip()
                if rec in record_types:
                    try:
                        x = float(line[30:38]) + translation[0]
                        y = float(line[38:46]) + translation[1]
                        z = float(line[46:54]) + translation[2]
                        out_lines.append(
                            line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
                        )
                    except ValueError:
                        out_lines.append(line)
                elif rec not in ("END", "TER", "CONECT"):
                    out_lines.append(line)

    out_lines = []
    _translate_pdb(protein_pdb, translation, ("ATOM",), out_lines)
    _translate_pdb(ligand_pdb,  translation, ("ATOM", "HETATM"), out_lines)
    out_lines.append("END\n")

    with open(output_path, "w") as f:
        f.writelines(out_lines)

    result["pdb_path"] = output_path
    result["centroid"] = centroid
    return result


def generate_tleap_script(protein_pdb:     str,
                           ligand_pdb:      str,
                           guest_prep:      str,
                           guest_frcmod:    str,
                           cofactor_params: list = None,
                           disulfide_pairs: list = None,
                           box_padding:     float = None,
                           output_top:      str   = "complex.top",
                           output_crd:      str   = "complex.inpcrd",
                           leap_script:     str   = "tleap_build.in") -> dict:
    """
    Generate tleap input script for protein-ligand system.

    Parameters
    ----------
    cofactor_params : list of {resname, prep_path, frcmod_path}
    disulfide_pairs : list of ((chainA, resA), (chainB, resB))
                      — tleap resids (sequential, from build_tleap_resid_map)
    box_padding     : float — solvation padding in Å (default 15.0)

    Returns
    -------
    dict with keys: script_path, warnings, errors
    """
    result = {
        "script_path": None,
        "warnings":    [],
        "errors":      [],
    }

    if box_padding is None:
        box_padding = SILENT_DEFAULTS["box_padding"]

    cofactor_params = cofactor_params or []
    disulfide_pairs = disulfide_pairs or []

    lines = [
        "# Auto-generated tleap script — LB-PaCS-MD workflow\n",
        "source leaprc.protein.ff19SB\n",
        "source leaprc.water.tip3p\n",
        "source leaprc.gaff2\n",
        "\n",
    ]

    # Load cofactor parameters
    for cof in cofactor_params:
        if cof.get("prep_path"):
            lines.append(f"loadAmberPrep {cof['prep_path']}\n")
        if cof.get("frcmod_path"):
            lines.append(f"loadAmberParams {cof['frcmod_path']}\n")

    # Load ligand parameters
    lines.append(f"loadAmberPrep {guest_prep}\n")
    lines.append(f"loadAmberParams {guest_frcmod}\n")
    lines.append("\n")

    # Load complex
    # Use combined centered PDB (protein + ligand together)
    complex_pdb = output_top.replace(".top", "_centered.pdb")
    if os.path.exists(ligand_pdb):
        # center_complex_at_origin was called before this
        complex_pdb = "complex_centered.pdb"

    lines.append(f"complex = loadpdb {complex_pdb}\n")
    lines.append("\n")

    # Disulfide bonds
    if disulfide_pairs:
        lines.append("# Disulfide bonds\n")
        for (ch_a, res_a), (ch_b, res_b) in disulfide_pairs:
            lines.append(
                f"bond complex.{res_a}.SG complex.{res_b}.SG\n"
            )
        lines.append("\n")

    # Solvation — octahedral box
    lines.append(
        f"solvateOct complex TIP3PBOX {box_padding:.1f}\n"
    )
    lines.append("\n")

    # Neutralize only
    lines.append("addions complex Na+ 0\n")
    lines.append("addions complex Cl- 0\n")
    lines.append("\n")

    # Save
    lines.append(f"saveamberparm complex {output_top} {output_crd}\n")
    lines.append("quit\n")

    with open(leap_script, "w") as f:
        f.writelines(lines)

    result["script_path"] = leap_script
    return result


def run_tleap(leap_script: str) -> dict:
    """
    Execute tleap with the given input script.

    Returns
    -------
    dict with keys:
        success       : bool
        log_path      : str
        n_atoms       : int  — system size (parsed from log)
        box_dims      : list — [a, b, c] in Å (parsed from log)
        warnings      : list
        errors        : list
    """
    result = {
        "success":   False,
        "log_path":  "tleap.log",
        "n_atoms":   None,
        "box_dims":  None,
        "warnings":  [],
        "errors":    [],
    }

    cmd = ["tleap", "-f", leap_script]

    try:
        with open("tleap.log", "w") as log:
            proc = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=300,
            )
    except FileNotFoundError:
        result["errors"].append(
            "tleap not found. Ensure AmberTools is installed."
        )
        return result
    except subprocess.TimeoutExpired:
        result["errors"].append("tleap timed out (>300s).")
        return result

    # Parse log for errors and system info
    with open("tleap.log", "r") as f:
        log_content = f.read()

    if "FATAL" in log_content or proc.returncode != 0:
        # Extract error lines
        err_lines = [l.strip() for l in log_content.splitlines()
                     if "FATAL" in l or "Error" in l or "error" in l]
        result["errors"].append(
            "tleap failed. Errors:\n" + "\n".join(err_lines[:10])
        )
        return result

    # Parse total atoms
    for line in log_content.splitlines():
        if "atoms" in line.lower() and "total" in line.lower():
            m = re.search(r"(\d+)\s+atoms", line)
            if m:
                result["n_atoms"] = int(m.group(1))

    # Parse box dimensions
    for line in log_content.splitlines():
        if "Box dimensions" in line or "box" in line.lower():
            nums = re.findall(r"[\d.]+", line)
            if len(nums) >= 3:
                result["box_dims"] = [float(n) for n in nums[:3]]
                break

    # Warn if box is very large
    if result["box_dims"] and max(result["box_dims"]) > 120:
        result["warnings"].append(
            f"Box dimensions are large ({result['box_dims']}). "
            "System may contain ~80-120k+ atoms. "
            "GPU memory requirements will be high."
        )

    result["success"] = True
    return result


# ============================================================
# SECTION J — MINIMIZATION + HEATING
# ============================================================

def run_minimization(top:         str,
                     inpcrd:      str,
                     temperature: float = 300.0,
                     output_rst7: str   = "last_frame.rst7") -> dict:
    """
    Quick minimization + heating using OpenMM.
    No restraints. Short protocol sufficient for BSV-placed ligand.

    Protocol (from SILENT_DEFAULTS):
        1. Minimize: max 5000 steps, tolerance 10 kJ/mol/nm
        2. Heat: 3 stages (T/3, 2T/3, T), 5000 steps each (10 ps each)
        3. NVT: 5000 steps (10 ps)
        4. Save: last_frame.rst7

    Parameters
    ----------
    top         : str   — AMBER topology (.top)
    inpcrd      : str   — AMBER coordinates (.inpcrd)
    temperature : float — target temperature in K
    output_rst7 : str   — output restart file

    Returns
    -------
    dict with keys:
        rst7_path    : str
        final_energy : float  — kJ/mol
        warnings     : list
        errors       : list
    """
    result = {
        "rst7_path":    None,
        "final_energy": None,
        "warnings":     [],
        "errors":       [],
    }

    try:
        import openmm as mm
        from openmm import app, unit, Platform
        import parmed as pmd
        import math as _math
    except ImportError as e:
        result["errors"].append(f"OpenMM/ParmEd import failed: {e}")
        return result

    # Detect best platform
    use_platform = _detect_openmm_platform()

    try:
        prmtop = app.AmberPrmtopFile(top)
        rst7   = app.AmberInpcrdFile(inpcrd)

        system = prmtop.createSystem(
            nonbondedMethod  = app.PME,
            nonbondedCutoff  = SILENT_DEFAULTS["cut"] * 0.1 * unit.nanometer,
            constraints      = app.HBonds,
            rigidWater       = True,
            ewaldErrorTolerance = SILENT_DEFAULTS["ewald_tolerance"],
        )

        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            SILENT_DEFAULTS["gamma_ln"] / unit.picosecond,
            SILENT_DEFAULTS["timestep_fs"] * unit.femtoseconds,
        )

        plat = Platform.getPlatformByName(use_platform)
        sim  = app.Simulation(prmtop.topology, system, integrator, plat)
        sim.context.setPositions(rst7.positions)

        if rst7.boxVectors is not None:
            sim.context.setPeriodicBoxVectors(*rst7.boxVectors)

        # Step 1: Minimize
        print(f"   Minimizing (max {SILENT_DEFAULTS['min_max_iterations']} steps)...")
        sim.minimizeEnergy(
            maxIterations = SILENT_DEFAULTS["min_max_iterations"],
            tolerance     = SILENT_DEFAULTS["min_tolerance"] * unit.kilojoule_per_mole / unit.nanometer,
        )

        # Step 2: Heat in stages
        stages = SILENT_DEFAULTS["heat_stages"]
        steps  = SILENT_DEFAULTS["heat_steps_per_stage"]
        temps  = [temperature * (i + 1) / stages for i in range(stages)]

        for T in temps:
            integrator.setTemperature(T * unit.kelvin)
            sim.context.setVelocitiesToTemperature(T * unit.kelvin)
            sim.step(steps)
            state = sim.context.getState(getEnergy=True)
            pe    = state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
            print(f"   Heated to {T:.0f} K: E = {pe:.1f} kJ/mol")

            if _math.isnan(pe) or _math.isinf(pe):
                result["errors"].append(
                    f"NaN/Inf energy detected during heating at {T:.0f} K."
                )
                return result

        # Step 3: NVT
        integrator.setTemperature(temperature * unit.kelvin)
        sim.step(SILENT_DEFAULTS["nvt_steps"])

        # Final energy
        state        = sim.context.getState(getEnergy=True, getPositions=True,
                                            getVelocities=True)
        final_energy = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)

        # Step 4: Save rst7
        save_rst7(sim, top, output_rst7)

        result["rst7_path"]    = output_rst7
        result["final_energy"] = round(final_energy, 2)

    except Exception as e:
        result["errors"].append(f"Minimization failed: {e}")

    return result


def _detect_openmm_platform() -> str:
    """
    Detect the best available OpenMM platform.
    Priority: OpenCL → CUDA → CPU.
    """
    try:
        from openmm import Platform
        platforms = [Platform.getPlatform(i).getName()
                     for i in range(Platform.getNumPlatforms())]
        for preferred in ("OpenCL", "CUDA"):
            if preferred in platforms:
                return preferred
    except Exception:
        pass
    return "CPU"


def save_rst7(simulation,
              top_path:    str,
              output_path: str = "last_frame.rst7") -> dict:
    """
    Save OpenMM simulation state as AMBER .rst7 restart file using ParmEd.

    Returns
    -------
    dict with keys: rst7_path, warnings, errors
    """
    result = {"rst7_path": None, "warnings": [], "errors": []}

    try:
        import parmed as pmd
        import numpy as np
        from openmm import unit
        import math

        state    = simulation.context.getState(
            getPositions=True, getVelocities=True
        )

        pos_nm   = state.getPositions(asNumpy=True)
        vel_nm_ps = state.getVelocities(asNumpy=True)

        pos_A    = pos_nm.value_in_unit(unit.nanometer) * 10.0
        vel_A_ps = vel_nm_ps.value_in_unit(unit.nanometer / unit.picosecond) * 10.0

        # Box vectors
        try:
            a_v, b_v, c_v = state.getPeriodicBoxVectors(asNumpy=True)
        except TypeError:
            a, b, c = state.getPeriodicBoxVectors()
            a_v = np.array([a.x, a.y, a.z])
            b_v = np.array([b.x, b.y, b.z])
            c_v = np.array([c.x, c.y, c.z])

        def vlen(v): return float(np.linalg.norm(v))
        def vang(u, v):
            cu = float(np.dot(u, v) / (vlen(u) * vlen(v)))
            cu = max(-1.0, min(1.0, cu))
            return math.degrees(math.acos(cu))

        A     = vlen(a_v) * 10.0
        B     = vlen(b_v) * 10.0
        C     = vlen(c_v) * 10.0
        alpha = vang(b_v, c_v)
        beta  = vang(a_v, c_v)
        gamma = vang(a_v, b_v)

        struct             = pmd.load_file(top_path)
        struct.coordinates = pos_A
        struct.velocities  = vel_A_ps
        struct.box         = [A, B, C, alpha, beta, gamma]
        struct.save(output_path, overwrite=True)

        result["rst7_path"] = output_path

    except Exception as e:
        result["errors"].append(f"save_rst7 failed: {e}")

    return result


# ============================================================
# SECTION K — PACS-MD OUTPUT + ZIP
# ============================================================

def generate_config_json(temperature:   float,
                          timestep:      float = None,
                          sim_time:      float = 20.0,
                          traj_interval: int   = 500,
                          output_path:   str   = "config.json") -> dict:
    """
    Generate OpenMM config.json for PaCS-Q.

    Returns
    -------
    dict with keys: config_path, warnings, errors
    """
    result = {"config_path": None, "warnings": [], "errors": []}

    if timestep is None:
        timestep = SILENT_DEFAULTS["timestep_fs"]

    time_step_ps  = timestep / 1000.0
    steps_cycle   = int(sim_time / time_step_ps)

    config = {
        "temperature":   str(temperature),
        "pressure":      str(SILENT_DEFAULTS["pres0"]),
        "timestep":      str(timestep),
        "friction":      str(SILENT_DEFAULTS["gamma_ln"]),
        "steps":         str(steps_cycle),
        "traj_interval": str(traj_interval),
    }

    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)

    result["config_path"] = output_path
    return result


def generate_mdin(temperature:   float,
                  timestep:      float = None,
                  sim_time:      float = 20.0,
                  traj_interval: int   = 500,
                  output_path:   str   = "md.mdin") -> dict:
    """
    Generate AMBER .mdin file for PaCS-Q.

    Returns
    -------
    dict with keys: mdin_path, warnings, errors
    """
    result = {"mdin_path": None, "warnings": [], "errors": []}

    if timestep is None:
        timestep = SILENT_DEFAULTS["timestep_fs"]

    dt_ps       = timestep / 1000.0
    nstlim      = int(sim_time / dt_ps)
    ntwx        = traj_interval

    mdin = f"""&cntrl
  imin   = 0,
  irest  = 1,
  ntx    = 5,
  nstlim = {nstlim},
  dt     = {dt_ps:.4f},
  ntc    = {SILENT_DEFAULTS['ntc']},
  ntf    = {SILENT_DEFAULTS['ntf']},
  cut    = {SILENT_DEFAULTS['cut']:.1f},
  ntt    = {SILENT_DEFAULTS['ntt']},
  gamma_ln = {SILENT_DEFAULTS['gamma_ln']:.1f},
  temp0  = {temperature:.1f},
  ntp    = {SILENT_DEFAULTS['ntp']},
  barostat = {SILENT_DEFAULTS['barostat']},
  pres0  = {SILENT_DEFAULTS['pres0']:.1f},
  ntwx   = {ntwx},
  ntwr   = {nstlim},
  ioutfm = 1,
  iwrap  = {SILENT_DEFAULTS['iwrap']},
&end
"""

    with open(output_path, "w") as f:
        f.write(mdin)

    result["mdin_path"] = output_path
    return result


def generate_mdp(temperature:   float,
                 timestep:      float = None,
                 sim_time:      float = 20.0,
                 traj_interval: int   = 500,
                 output_path:   str   = "md.mdp") -> dict:
    """
    Generate GROMACS .mdp file.
    NOTE: PaCS-Q does NOT support GROMACS — this is for manual workflow only.

    Returns
    -------
    dict with keys: mdp_path, warnings, errors
    """
    result = {"mdp_path": None, "warnings": [], "errors": []}
    result["warnings"].append(
        "GROMACS output is for manual workflow only. "
        "PaCS-Q does not support GROMACS engine."
    )

    if timestep is None:
        timestep = SILENT_DEFAULTS["timestep_fs"]

    dt_ps  = timestep / 1000.0
    nsteps = int(sim_time / dt_ps)

    mdp = f"""integrator  = md
nsteps      = {nsteps}
dt          = {dt_ps:.4f}
nstxout     = {traj_interval}
nstvout     = {traj_interval}
nstfout     = 0
nstlog      = 500
nstxout-compressed = {traj_interval}

; Bond parameters
continuation        = yes
constraint-algorithm = lincs
constraints         = h-bonds
lincs-iter          = 1
lincs-order         = 4

; Nonbonded settings
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres

; Electrostatics
coulombtype     = PME
pme-order       = 4
fourierspacing  = 0.16

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = {temperature:.1f} {temperature:.1f}

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = {SILENT_DEFAULTS['pres0']:.1f}
compressibility = 4.5e-5

; Periodic boundary conditions
pbc = xyz
"""

    with open(output_path, "w") as f:
        f.write(mdp)

    result["mdp_path"] = output_path
    return result


def amber_to_gromacs(top:           str,
                     inpcrd:        str,
                     rst7:          str,
                     output_prefix: str = "complex_gmx") -> dict:
    """
    Convert AMBER topology/coordinates to GROMACS format using ParmEd.
    Falls back to acpype if ParmEd fails.

    Returns
    -------
    dict with keys:
        gro_path  : str
        top_path  : str
        itp_path  : str (ligand ITP if extracted)
        warnings  : list
        errors    : list
    """
    result = {
        "gro_path":  None,
        "top_path":  None,
        "itp_path":  None,
        "warnings":  [],
        "errors":    [],
    }
    result["warnings"].append(
        "GROMACS conversion is for reference only. "
        "PaCS-Q does not support GROMACS."
    )

    try:
        import parmed as pmd
        struct = pmd.load_file(top, rst7)

        gro_path = f"{output_prefix}.gro"
        top_path = f"{output_prefix}.top"

        struct.save(gro_path, overwrite=True)
        struct.save(top_path, overwrite=True)

        result["gro_path"] = gro_path
        result["top_path"] = top_path

    except Exception as e:
        result["errors"].append(f"ParmEd conversion failed: {e}")

        # Fallback: try acpype
        try:
            cmd = ["acpype", "-p", top, "-x", inpcrd, "-b", output_prefix]
            subprocess.run(cmd, check=True, timeout=120)
            result["warnings"].append("Used acpype fallback for GROMACS conversion.")
        except Exception as e2:
            result["errors"].append(f"acpype fallback also failed: {e2}")

    return result


def validate_pacs_params(cycles:     int,
                          candidates: int,
                          sim_time:   float) -> dict:
    """
    Validate PaCS-MD parameters and return warnings.

    Returns
    -------
    dict with keys: valid, warnings, errors
    """
    result = {"valid": True, "warnings": [], "errors": []}

    if cycles < 10:
        result["warnings"].append(
            f"cycles={cycles} is low. Ligand may not reach the binding site. "
            "Recommended: ≥ 10 cycles."
        )
    if candidates < 3:
        result["warnings"].append(
            f"candidates={candidates} may give poor pathway sampling. "
            "Recommended: ≥ 3."
        )
    if sim_time < 10:
        result["warnings"].append(
            f"sim_time={sim_time} ps/cycle is very short. "
            "Recommended: ≥ 10 ps/cycle."
        )
    total_ps = cycles * candidates * sim_time
    if total_ps > 50000:
        result["warnings"].append(
            f"Total simulation ≈ {total_ps/1000:.1f} ns. "
            "Consider running on HPC."
        )

    return result


def estimate_compute_time(cycles:     int,
                           candidates: int,
                           sim_time:   float,
                           engine:     str,
                           n_atoms:    int = 80000) -> dict:
    """
    Rough estimate of compute time for LB-PaCS-MD.

    Based on approximate benchmarks:
        CPU    : ~0.5 ns/day per 10k atoms
        GPU T4 : ~5 ns/day per 10k atoms
        GPU A100: ~20 ns/day per 10k atoms

    Returns
    -------
    dict with keys: cpu_hours, gpu_t4_minutes, gpu_a100_minutes
    """
    total_ns   = cycles * candidates * sim_time / 1000.0
    atom_scale = n_atoms / 10000.0

    cpu_ns_per_day   = 0.5 / atom_scale
    t4_ns_per_day    = 5.0 / atom_scale
    a100_ns_per_day  = 20.0 / atom_scale

    return {
        "total_ns":        round(total_ns, 2),
        "cpu_hours":       round(total_ns / cpu_ns_per_day * 24, 1),
        "gpu_t4_minutes":  round(total_ns / t4_ns_per_day * 24 * 60, 1),
        "gpu_a100_minutes":round(total_ns / a100_ns_per_day * 24 * 60, 1),
    }


def generate_run_script(engine:          str,
                         top:             str,
                         rst7:            str,
                         cycles:          int,
                         candidates:      int,
                         binding_site_sel: str,
                         temperature:     float = 300.0,
                         sim_time:        float = 20.0,
                         traj_interval:   int   = 500,
                         scheduler:       str   = None,
                         scheduler_params: dict  = None,
                         output_path:     str   = "run_pacs_md.sh") -> dict:
    """
    Generate run_pacs_md.sh shell script.

    Parameters
    ----------
    engine          : str  — "openmm" | "amber"
    scheduler       : str  — None | "slurm" | "pbs" | "sge"
    scheduler_params: dict — {partition, nodes, cpus, walltime, module}

    Returns
    -------
    dict with keys: script_path, warnings, errors
    """
    result = {"script_path": None, "warnings": [], "errors": []}

    sp = scheduler_params or {}

    # Build pacs_q_md command
    engine_flag = engine.lower()
    if engine_flag == "amber":
        engine_args = f"-e amber"
    else:
        engine_args = f"-e openmm -cf config.json"

    pacs_cmd = (
        f"pacs_q_md \\\n"
        f"  -p {top} \\\n"
        f"  -c {rst7} \\\n"
        f"  -cy {cycles} \\\n"
        f"  -cd {candidates} \\\n"
        f"  -s \"{binding_site_sel}\" \\\n"
        f"  -s2 \"resname MOL\" \\\n"
        f"  -b \\\n"
        f"  {engine_args}"
    )

    # Local script
    local_lines = [
        "#!/bin/bash\n",
        "# LB-PaCS-MD run script — auto-generated by LB-PaCS-MD workflow\n",
        f"# Engine: {engine} | Cycles: {cycles} | Candidates: {candidates}\n",
        f"# Temperature: {temperature} K | Sim/cycle: {sim_time} ps\n",
        "\n",
        "# Check / install PaCS-Q\n",
        "if ! command -v pacs_q_md &>/dev/null; then\n",
        "    echo 'PaCS-Q not found — installing...'\n",
        "    pip install PaCS-Q\n",
        "fi\n",
        "\n",
        pacs_cmd + "\n",
    ]

    if scheduler == "slurm":
        partition = sp.get("partition", "gpu")
        nodes     = sp.get("nodes", 1)
        cpus      = sp.get("cpus", 4)
        walltime  = sp.get("walltime", "24:00:00")
        module    = sp.get("module", "")

        slurm_lines = [
            "#!/bin/bash\n",
            f"#SBATCH --job-name=LB_PaCS_MD\n",
            f"#SBATCH --partition={partition}\n",
            f"#SBATCH --nodes={nodes}\n",
            f"#SBATCH --cpus-per-task={cpus}\n",
            f"#SBATCH --gres=gpu:1\n",
            f"#SBATCH --time={walltime}\n",
            f"#SBATCH --output=pacs_md_%j.log\n",
            "\n",
        ]
        if module:
            slurm_lines += [f"module load {module}\n", "\n"]
        slurm_lines += [
            "cd $SLURM_SUBMIT_DIR\n",
            "\n",
            "# Check / install PaCS-Q\n",
            "if ! command -v pacs_q_md &>/dev/null; then\n",
            "    pip install PaCS-Q\n",
            "fi\n",
            "\n",
            pacs_cmd + "\n",
        ]
        script_lines = slurm_lines

    elif scheduler == "pbs":
        walltime = sp.get("walltime", "24:00:00")
        module   = sp.get("module", "")

        pbs_lines = [
            "#!/bin/bash\n",
            f"#PBS -N LB_PaCS_MD\n",
            f"#PBS -l walltime={walltime}\n",
            f"#PBS -l select=1:ngpus=1\n",
            f"#PBS -o pacs_md.log\n",
            "\n",
        ]
        if module:
            pbs_lines += [f"module load {module}\n", "\n"]
        pbs_lines += [
            "cd $PBS_O_WORKDIR\n",
            "\n",
            "if ! command -v pacs_q_md &>/dev/null; then\n",
            "    pip install PaCS-Q\n",
            "fi\n",
            "\n",
            pacs_cmd + "\n",
        ]
        script_lines = pbs_lines

    else:
        script_lines = local_lines

    with open(output_path, "w") as f:
        f.writelines(script_lines)
    os.chmod(output_path, 0o755)

    result["script_path"] = output_path
    return result


def generate_summary(all_params: dict,
                     output_path: str = "summary.txt") -> dict:
    """
    Generate human-readable summary of the entire workflow.

    Parameters
    ----------
    all_params : dict — collected from all cells (see global variables flow)

    Returns
    -------
    dict with keys: summary_path, warnings, errors
    """
    result = {"summary_path": None, "warnings": [], "errors": []}

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 60 + "\n",
        "LB-PaCS-MD Workflow Summary\n",
        f"Generated: {now}\n",
        "=" * 60 + "\n",
        "\n",
        "[ RECEPTOR ]\n",
        f"  PDB source         : {all_params.get('pdb_source', 'N/A')}\n",
        f"  Protein pH         : {all_params.get('target_pH', 'N/A')}\n",
        f"  Binding site resids: {all_params.get('binding_site_sel', 'N/A')}\n",
        f"  Disulfide pairs    : {all_params.get('disulfide_pairs', [])}\n",
        "\n",
        "[ LIGAND ]\n",
        f"  Input mode         : {all_params.get('ligand_mode', 'N/A')}\n",
        f"  Residue name       : MOL\n",
        f"  Charge             : {all_params.get('ligand_charge', 'N/A')}\n",
        "\n",
        "[ BSV PLACEMENT ]\n",
        f"  Placement distance : {all_params.get('placement_distance', 'N/A')} Å\n",
        f"  Buffer             : {all_params.get('buffer', SILENT_DEFAULTS['default_buffer'])} Å\n",
        f"  Exit direction     : {all_params.get('exit_direction', 'N/A')}\n",
        "\n",
        "[ TOPOLOGY ]\n",
        f"  Force field        : ff19SB + GAFF2 + TIP3P\n",
        f"  Box type           : Octahedral\n",
        f"  Box padding        : {SILENT_DEFAULTS['box_padding']} Å\n",
        f"  System atoms       : {all_params.get('n_atoms', 'N/A')}\n",
        "\n",
        "[ MINIMIZATION ]\n",
        f"  Temperature        : {all_params.get('target_temperature', 300)} K\n",
        f"  Protocol           : Minimize → Heat (3 stages) → NVT 10ps\n",
        "\n",
        "[ LB-PaCS-MD ]\n",
        f"  Engine             : {all_params.get('engine', 'N/A')}\n",
        f"  Cycles             : {all_params.get('cycles', 'N/A')}\n",
        f"  Candidates/cycle   : {all_params.get('candidates', 'N/A')}\n",
        f"  Sim time/cycle     : {all_params.get('sim_time', 'N/A')} ps\n",
        f"  Timestep           : {SILENT_DEFAULTS['timestep_fs']} fs\n",
        f"  Mode               : Binding (-b)\n",
        "\n",
        "[ FILES ]\n",
        f"  Topology           : system/complex.top\n",
        f"  Coordinates        : system/complex.inpcrd\n",
        f"  Restart            : system/last_frame.rst7\n",
        f"  Run script         : run_pacs_md.sh\n",
        "\n",
        "=" * 60 + "\n",
        "Reference: Kuvek et al., J. Chem. Inf. Model. 2026, 66, 1604-1616\n",
        "           DOI: 10.1021/acs.jcim.5c02705\n",
        "=" * 60 + "\n",
    ]

    with open(output_path, "w") as f:
        f.writelines(lines)

    result["summary_path"] = output_path
    return result


def generate_readme(system_info: dict,
                    engine:      str,
                    scheduler:   str  = None,
                    output_path: str  = "README.md") -> dict:
    """
    Generate README.md for the downloaded zip package.

    Returns
    -------
    dict with keys: readme_path, warnings, errors
    """
    result = {"readme_path": None, "warnings": [], "errors": []}

    engine_note = {
        "openmm": "Requires OpenMM. Install: `conda install -c conda-forge openmm`",
        "amber":  "Requires AMBER ≥ 22. Ensure `pmemd` or `sander` is in PATH.",
    }.get(engine.lower(), "")

    sched_note = ""
    if scheduler == "slurm":
        sched_note = "\nSubmit to SLURM: `sbatch run_pacs_md.sh`"
    elif scheduler == "pbs":
        sched_note = "\nSubmit to PBS: `qsub run_pacs_md.sh`"
    else:
        sched_note = "\nRun locally: `bash run_pacs_md.sh`"

    content = f"""# LB-PaCS-MD Input Package

Generated by LB-PaCS-MD workflow — {datetime.now().strftime("%Y-%m-%d")}

## Contents

```
system/
├── complex.top       AMBER topology
├── complex.inpcrd    Initial coordinates
└── last_frame.rst7   Minimized/heated restart

components/
├── protein_prepared.pdb
├── ligand_placed.pdb
├── guest.prep
└── guest.frcmod

config.json / md.mdin   Engine configuration
run_pacs_md.sh          Run script
summary.txt             Workflow summary
README.md               This file
```

## Prerequisites

- PaCS-Q: `pip install PaCS-Q`
- {engine_note}

## Running LB-PaCS-MD
{sched_note}

## Post-Analysis

Post-analysis (trajectory visualization, binding pathway analysis,
binding free energy) is the user's responsibility.

Recommended tools:
- MDAnalysis (Python)
- CPPTRAJ (AmberTools)
- VMD

## Reference

Kuvek et al., J. Chem. Inf. Model. 2026, 66, 1604-1616
DOI: 10.1021/acs.jcim.5c02705
"""

    with open(output_path, "w") as f:
        f.write(content)

    result["readme_path"] = output_path
    return result


def create_zip(engine:          str,
               top_path:        str,
               inpcrd_path:     str,
               rst7_path:       str,
               protein_pdb:     str,
               ligand_pdb:      str,
               guest_prep:      str,
               guest_frcmod:    str,
               run_script:      str,
               summary_path:    str,
               readme_path:     str,
               config_or_mdin:  str = None,
               output_zip:      str = "LB_PaCS_MD_inputs.zip") -> dict:
    """
    Create download zip package with all LB-PaCS-MD input files.

    Zip structure:
        system/
            complex.top
            complex.inpcrd
            last_frame.rst7
        components/
            protein_prepared.pdb
            ligand_placed.pdb
            guest.prep
            guest.frcmod
        config.json / md.mdin
        run_pacs_md.sh
        summary.txt
        README.md

    Returns
    -------
    dict with keys: zip_path, warnings, errors
    """
    result = {"zip_path": None, "warnings": [], "errors": []}

    files_map = {
        "system/complex.top":           top_path,
        "system/complex.inpcrd":        inpcrd_path,
        "system/last_frame.rst7":       rst7_path,
        "components/protein_prepared.pdb": protein_pdb,
        "components/ligand_placed.pdb": ligand_pdb,
        "components/guest.prep":        guest_prep,
        "components/guest.frcmod":      guest_frcmod,
        "run_pacs_md.sh":               run_script,
        "summary.txt":                  summary_path,
        "README.md":                    readme_path,
    }

    if config_or_mdin:
        basename = os.path.basename(config_or_mdin)
        files_map[basename] = config_or_mdin

    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for arcname, filepath in files_map.items():
            if filepath and os.path.exists(filepath):
                zf.write(filepath, arcname)
            else:
                result["warnings"].append(
                    f"File not found, skipped in zip: {filepath} ({arcname})"
                )

    result["zip_path"] = output_zip
    return result


def validate_inputs_for_pacs(top_path:     str,
                               rst7_path:    str,
                               engine:       str,
                               binding_sel:  str,
                               temperature:  float,
                               cycles:       int,
                               candidates:   int,
                               sim_time:     float) -> dict:
    """
    Final validation before generating PaCS-MD inputs.
    Performs all checks and returns traffic-light style report.

    Returns
    -------
    dict with keys:
        passed    : bool   — True only if no errors
        checks    : list of {name, status ("ok"|"warn"|"error"), message}
        warnings  : list
        errors    : list
    """
    result = {
        "passed":   True,
        "checks":   [],
        "warnings": [],
        "errors":   [],
    }

    def _add(name, status, msg):
        result["checks"].append({
            "name":    name,
            "status":  status,
            "message": msg,
        })
        if status == "error":
            result["errors"].append(msg)
            result["passed"] = False
        elif status == "warn":
            result["warnings"].append(msg)

    # File checks
    _add("complex.top",
         "ok" if os.path.exists(top_path) else "error",
         f"Topology: {top_path}" if os.path.exists(top_path)
         else f"Topology NOT FOUND: {top_path}")

    _add("last_frame.rst7",
         "ok" if os.path.exists(rst7_path) else "error",
         f"Restart: {rst7_path}" if os.path.exists(rst7_path)
         else f"Restart NOT FOUND: {rst7_path}")

    # PaCS-Q check
    which = subprocess.run(["which", "pacs_q_md"],
                           capture_output=True, text=True)
    _add("PaCS-Q",
         "ok" if which.returncode == 0 else "warn",
         "PaCS-Q found: " + which.stdout.strip()
         if which.returncode == 0
         else "PaCS-Q not found — will be installed via pip at runtime")

    # Engine check
    if engine.lower() == "openmm":
        try:
            import openmm
            _add("OpenMM", "ok",
                 f"OpenMM {openmm.__version__} available")
        except ImportError:
            _add("OpenMM", "error",
                 "OpenMM not installed. Run: conda install -c conda-forge openmm")
    elif engine.lower() == "amber":
        amber_check = subprocess.run(
            ["which", "sander"], capture_output=True, text=True)
        _add("AMBER",
             "ok" if amber_check.returncode == 0 else "error",
             "sander found" if amber_check.returncode == 0
             else "sander not found. Ensure AMBER ≥ 22 is installed and in PATH.")

    # Selection check
    _add("Protein selection",
         "ok" if binding_sel and len(binding_sel) > 0 else "error",
         f"Protein selection: {binding_sel}" if binding_sel
         else "Protein selection is empty")

    _add("Ligand selection", "ok", "resname MOL (auto)")

    # Temperature from Cell 14
    _add("Temperature", "ok",
         f"Temperature: {temperature} K (from minimization cell)")

    # Parameter warnings
    param_check = validate_pacs_params(cycles, candidates, sim_time)
    for w in param_check["warnings"]:
        _add("Parameters", "warn", w)
    if not param_check["warnings"]:
        _add("Parameters", "ok",
             f"cycles={cycles}, candidates={candidates}, sim_time={sim_time} ps")

    return result
