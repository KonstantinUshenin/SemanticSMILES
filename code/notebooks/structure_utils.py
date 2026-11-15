"""
structure_utils.py
------------------

Helper utilities for SMILES structure classification.
Used by compute_Structures() in cheminformatics datasets.

Functions provided:
- split_fragments: safely split RDKit Mol into fragments
- get_fragment_charge: return formal charge of a fragment
- is_transition_metal: detect single-atom transition metal
- is_protonated_positive / is_protonated_negative: detect protonated charged species
- is_zwitterion: detect internal + and - within one molecule
- detect_reaction_smiles: detect reaction-style SMILES with ">>" or ">"
- detect_solvents: identify common solvents or hydrates
"""

from typing import List
from rdkit import Chem
from rdkit.Chem import rdchem
from typing import Set, Dict, Tuple


pt = rdchem.GetPeriodicTable()

def split_fragments(mol: Chem.Mol) -> List[Chem.Mol]:
    """Split RDKit molecule into disconnected fragments (ions, solvates, etc.)."""
    if mol is None:
        return []
    return Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

def get_fragment_charge(f: Chem.Mol) -> int:
    """Return formal charge of a fragment."""
    try:
        return Chem.GetFormalCharge(f)
    except Exception:
        return 0

def is_transition_metal(f: Chem.Mol) -> bool:
    """Return True if fragment is a single-atom transition metal (neutral or charged)."""
    metals = [pt.GetElementSymbol(i) for i in range(21, 31)]  # Sc–Zn
    metals += [pt.GetElementSymbol(i) for i in range(39, 49)]  # Y–Cd
    metals += [pt.GetElementSymbol(i) for i in range(72, 81)]  # Hf–Hg
    atoms = [a.GetSymbol() for a in f.GetAtoms()]
    return len(atoms) == 1 and atoms[0] in metals

def is_organic(f: Chem.Mol) -> bool:
    """Return True if fragment contains carbon and no transition metals."""
    symbols = [a.GetSymbol() for a in f.GetAtoms()]
    if "C" not in symbols:
        return False
    return not any(is_transition_metal(Chem.MolFromSmiles(a.GetSymbol())) for a in f.GetAtoms())

def is_protonated_positive(f: Chem.Mol) -> bool:
    """Return True if molecule is positively charged and contains hydrogen (e.g., protonated amine)."""
    atoms = [a.GetSymbol() for a in f.GetAtoms()]
    q = Chem.GetFormalCharge(f)
    if atoms == ["H"] and q == 1:
        return True
    return q > 0 and "H" in atoms

def is_protonated_negative(f: Chem.Mol) -> bool:
    """Return True if molecule is negatively charged and contains hydrogen (e.g., hydride, deprotonated acid)."""
    atoms = [a.GetSymbol() for a in f.GetAtoms()]
    q = Chem.GetFormalCharge(f)
    if atoms == ["H"] and q == -1:
        return True
    return q < 0 and "H" in atoms

def is_zwitterion(f: Chem.Mol) -> bool:
    """Return True if molecule contains both positive and negative atoms but total charge = 0."""
    has_pos = any(a.GetFormalCharge() > 0 for a in f.GetAtoms())
    has_neg = any(a.GetFormalCharge() < 0 for a in f.GetAtoms())
    return has_pos and has_neg and Chem.GetFormalCharge(f) == 0

def detect_reaction_smiles(smiles: str) -> bool:
    """Detect if SMILES represents a reaction or equilibrium (contains '>' or '>>')."""
    return ">>" in smiles or ">" in smiles

COMMON_SOLVENTS = {
    "O", "CO", "CCO", "CCCO", "C(C)O", "C=O", "CC(=O)O", "C(C)OC", "CCOCC", "CCCOC", "OCCO"
}

def detect_solvents(fragments: List[Chem.Mol]) -> int:
    """Count how many fragments correspond to common solvents."""
    from rdkit import Chem
    count = 0
    for f in fragments:
        s = Chem.MolToSmiles(f)
        if s in COMMON_SOLVENTS:
            count += 1
    return count

TOOLS_RULES = {
    # --- ANI family ---
    "ANI-1": {
        "allowed_elements": {"H", "C", "N", "O"},
        "allow_charged": False,
        "description": "ANI-1/1x/1ccx — high-accuracy neural potentials trained on DFT and CCSD(T) data for neutral organic molecules (H, C, N, O).",
    },
    "ANI-2x": {
        "allowed_elements": {"H", "C", "N", "O", "F", "Cl", "S"},
        "allow_charged": False,
        "description": "ANI-2x/3x — extended ANI model covering halogens and sulfur; neutral drug-like molecules.",
    },
    "Schrodinger-NeuroChem-ANI": {
        "allowed_elements": {"H", "C", "N", "O", "S", "F", "Cl", "P"},
        "allow_charged": False,
        "description": "Schrödinger-NeuroChem ANI — proprietary network for H, C, N, O, S, F, Cl, P; neutral organic and drug-like chemistry.",
    },

    # --- MACE family ---
    "MACE-ANI-CC": {
        "allowed_elements": {"H", "C", "N", "O"},
        "allow_charged": False,
        "description": "MACE-ANI-CC — coupled-cluster-accurate potential for neutral organics (H, C, N, O).",
    },
    "MACE-OFF23": {
        "allowed_elements": {"H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"},
        "allow_charged": False,
        "description": "MACE-OFF23 — transferable organic force field (10 elements) for neutral molecules.",
    },

    # --- Neural DFT surrogate ---
    "nabla2DFT": {
        "allowed_elements": {"H", "C", "N", "S", "O", "F", "Cl", "Br"},
        "allow_charged": False,
        "description": "nabla2DFT — neural DFT model for neutral organic systems; trained on H, C, N, O, S, F, Cl, Br.",
    },

    # --- Semiempirical / hybrid MLIPs ---
    "xTB": {
        "max_atomic_number": 86,
        "allow_charged": True,
        "description": "GFN-xTB family — semiempirical tight-binding for elements up to Z ≤ 86; supports charged and neutral systems.",
    },
    "AIQM1": {
        "allowed_elements": {"H", "C", "N", "O"},
        "allow_charged": True,
        "description": "AIQM1 — hybrid xTB + neural correction; DFT-level accuracy for charged and neutral organics (H, C, N, O).",
    },
    "AIMNet": {
        "allowed_elements": {"H", "C", "N", "O", "F", "Cl", "S"},
        "allow_charged": True,
        "description": "AIMNet / AIMNet2 / AIMNet3 — multi-task NNP predicting energy, charge, and dipoles; supports charged molecules.",
    },
    "FENNIX": {
        "allowed_elements": {"H", "C", "N", "O", "F", "S"},
        "allow_charged": True,
        "description": "FENNIX — hybrid ML + xTB correction model for charged and reactive organic chemistry.",
    },
}


# -------------------------------------------------------------
# Basic extractors
# -------------------------------------------------------------

def get_elements(mol: Chem.Mol) -> Set[str]:
    """Return a set of element symbols in the molecule."""
    return {a.GetSymbol() for a in mol.GetAtoms()} if mol else set()

def get_atomic_numbers(mol: Chem.Mol) -> Set[int]:
    """Return a set of atomic numbers (Z) in the molecule."""
    return {int(a.GetAtomicNum()) for a in mol.GetAtoms()} if mol else set()

def get_total_charge(mol: Chem.Mol) -> int:
    """Return total formal charge of the molecule."""
    try:
        return Chem.GetFormalCharge(mol)
    except Exception:
        return 0

def is_neutral(mol: Chem.Mol) -> bool:
    """Return True if molecule is neutral."""
    return get_total_charge(mol) == 0

# -------------------------------------------------------------
# Predicates
# -------------------------------------------------------------

def elements_within(elements: Set[str], allowed: Set[str]) -> bool:
    """Return True if all elements are within the allowed set."""
    return all(e in allowed for e in elements)

def atomic_numbers_within_max(atomic_numbers: Set[int], max_Z: int) -> bool:
    """Return True if all atomic numbers are <= max_Z."""
    return all(1 <= z <= max_Z for z in atomic_numbers)

# -------------------------------------------------------------
# Tool evaluation
# -------------------------------------------------------------

def check_tool_readiness(mol: Chem.Mol, tool: str) -> Tuple[bool, str]:
    """
    Check readiness for a single tool.
    Returns (ready_bool, reason_if_not_ready_or_empty_string).
    """
    rules = TOOLS_RULES[tool]
    charge = get_total_charge(mol)

    # xTB: Z <= 86, charge allowed
    if tool == "xTB":
        Zs = get_atomic_numbers(mol)
        if not Zs:
            return False, "empty_molecule"
        if not atomic_numbers_within_max(Zs, int(rules["max_atomic_number"])):
            return False, f"contains_Z_gt_{rules['max_atomic_number']}"
        return True, ""

    # List-based tools: must be neutral + element subset
    elems = get_elements(mol)
    if not rules.get("allow_charged", False) and charge != 0:
        return False, "non_neutral_molecule"
    allowed = rules["allowed_elements"]
    if not elements_within(elems, allowed):
        missing = sorted(e for e in elems if e not in allowed)
        return False, "unsupported_elements:" + ",".join(missing)
    return True, ""