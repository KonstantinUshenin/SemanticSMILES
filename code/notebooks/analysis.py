from typing import Dict, Any, List, Callable, OrderedDict, Optional, Union, Set
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    # Package-style import (e.g., `import notebooks.analysis`)
    from .structure_utils import (
        split_fragments,
        get_fragment_charge,
        is_transition_metal,
        is_protonated_positive,
        is_protonated_negative,
        is_zwitterion,
        detect_reaction_smiles,
        detect_solvents,
        is_organic
    )
except ImportError:
    # Notebook-folder import (e.g., cwd=notebooks, `import analysis`)
    from structure_utils import (
        split_fragments,
        get_fragment_charge,
        is_transition_metal,
        is_protonated_positive,
        is_protonated_negative,
        is_zwitterion,
        detect_reaction_smiles,
        detect_solvents,
        is_organic
    )


# --- Global registry for parameter functions ---
PARAMETER_REGISTRY: Dict[str, Dict[str, Any]] = OrderedDict()


def compute_parameter(name: str,
                      description: str,
                      aggregation: List[str] = []):
    """
    Decorator that registers a molecular property computation function
    with metadata.
    """
    def decorator(func: Callable[[str, Any, Dict[str, Any], Dict[str, Any]], Any]):
        PARAMETER_REGISTRY[name] = {
            "name": name,
            "description": description,
            "aggregation": aggregation,
            "function": func
        }
        return func
    return decorator


# @compute_parameter(
#     "molecular_weight",
#     "Molecular weight",
#     ["_mean", "_min", "_max"],
#     ["sum", "min", "max"],
#     ["average", None, None]
# )
# def compute_MolWt(smiles: str,
#                   mol: Any,
#                   row_raw: Dict[str, object],
#                   row_processed: Dict[str, object]):
#     """Compute molecular weight from RDKit Mol."""
#     if mol is None:
#         return None
#     try:
#         return Descriptors.MolWt(mol)
#     except Exception:
#         return None

@compute_parameter("molecular_weight",
                   "Molecular weight (Daltons)",
                   ['min','mean','max'])
def compute_MolWt(smiles: str, mol: Any,
                  row_raw: Dict[str, object],
                  row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.MolWt(mol)
    except Exception:
        return None


@compute_parameter("mol_logp", "Octanol–water partition coefficient (logP)")
def compute_MolLogP(smiles: str, mol: Any,
                    row_raw: Dict[str, object],
                    row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.MolLogP(mol)
    except Exception:
        return None


@compute_parameter("mol_mr", "Molar refractivity")
def compute_MolMR(smiles: str, mol: Any,
                  row_raw: Dict[str, object],
                  row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.MolMR(mol)
    except Exception:
        return None


@compute_parameter("num_heteroatoms", "Number of heteroatoms (non-C/H)")
def compute_NumHeteroatoms(smiles: str, mol: Any,
                           row_raw: Dict[str, object],
                           row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.NumHeteroatoms(mol)
    except Exception:
        return None


@compute_parameter("heavy_atom_count", "Number of heavy (non-H) atoms")
def compute_HeavyAtomCount(smiles: str, mol: Any,
                           row_raw: Dict[str, object],
                           row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.HeavyAtomCount(mol)
    except Exception:
        return None


@compute_parameter("ring_count", "Number of rings in molecule")
def compute_RingCount(smiles: str, mol: Any,
                      row_raw: Dict[str, object],
                      row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.RingCount(mol)
    except Exception:
        return None


@compute_parameter("tpsa", "Topological Polar Surface Area (Å²)")
def compute_TPSA(smiles: str, mol: Any,
                 row_raw: Dict[str, object],
                 row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.TPSA(mol)
    except Exception:
        return None


@compute_parameter("num_h_donors", "Number of hydrogen bond donors")
def compute_NumHDonors(smiles: str, mol: Any,
                       row_raw: Dict[str, object],
                       row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.NumHDonors(mol)
    except Exception:
        return None


@compute_parameter("num_h_acceptors", "Number of hydrogen bond acceptors")
def compute_NumHAcceptors(smiles: str, mol: Any,
                          row_raw: Dict[str, object],
                          row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.NumHAcceptors(mol)
    except Exception:
        return None


@compute_parameter("fraction_csp3", "Fraction of sp³ hybridized carbons")
def compute_FractionCSP3(smiles: str, mol: Any,
                         row_raw: Dict[str, object],
                         row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.FractionCSP3(mol)
    except Exception:
        return None


@compute_parameter("num_rotatable_bonds", "Number of rotatable bonds")
def compute_NumRotatableBonds(smiles: str, mol: Any,
                              row_raw: Dict[str, object],
                              row_processed: Dict[str, object]):
    if mol is None:
        return None
    try:
        return Descriptors.NumRotatableBonds(mol)
    except Exception:
        return None

@compute_parameter("elements", "Set of elements in molecule", ["unionset", "uniondict"])
def compute_Elements(smiles: str,
                     mol: Any,
                     row_raw: Dict[str, object],
                     row_processed: Dict[str, object]):
    """Return the set of unique elements present in the molecule."""
    if mol is None:
        return None
    try:
        return {atom.GetSymbol() for atom in mol.GetAtoms()}
    except Exception:
        return None


@compute_parameter("elementsinmol", "Element count dictionary", ["uniondict"])
def compute_Elements(smiles: str,
                     mol: Any,
                     row_raw: Dict[str, object],
                     row_processed: Dict[str, object]):
    """Return a dictionary {element_symbol: count} for the given molecule."""
    if mol is None:
        return None
    try:
        element_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in element_counts:
                element_counts[symbol] += 1
            else:
                element_counts[symbol] = 1
        return element_counts
    except Exception:
        return None


@compute_parameter("structures", "Type of structures in records", ["uniondict"])
def compute_Structures(smiles: str,
                     mol: Any,
                     row_raw: Dict[str, object],
                     row_processed: Dict[str, object]):
    """Return a dictionary {structure_type: count} for the given molecule."""
    if mol is None:
        return None
    try:
        structures_counts = {}
        # The code here
        if False:
            pass
        else:
            structures_counts = {"other": 1}
        return structures_counts
    except Exception:
        return None

@compute_parameter("structures", "Type of structures in records", ["uniondict"])
def compute_Structures(smiles: str,
                       mol: Any,
                       row_raw: Dict[str, object],
                       row_processed: Dict[str, object]):
    """Classify molecule(s) by structural type.
    Returns {structure_type: 1} dictionary.
    """

    if mol is None:
        return None

    try:
        # Reaction SMILES
        if detect_reaction_smiles(smiles):
            return {"reaction_smiles": 1}

        fragments = split_fragments(mol)
        n_frags = len(fragments)
        frag_charges = [get_fragment_charge(f) for f in fragments]
        total_charge = sum(frag_charges)
        structures_counts = {}

        # ---- Single fragment ----
        if n_frags == 1:
            f = fragments[0]
            q = frag_charges[0]

            if is_transition_metal(f):
                structures_counts["transition_metal_neutral" if q == 0 else "transition_metal_complex"] = 1
            elif is_zwitterion(f):
                structures_counts["zwitterion"] = 1
            elif q == 0:
                structures_counts["neutral_molecule"] = 1
            elif q > 0 and is_protonated_positive(f):
                structures_counts["protonated_positive"] = 1
            elif q < 0 and is_protonated_negative(f):
                structures_counts["protonated_negative"] = 1
            elif q > 0:
                structures_counts["organic_cation" if is_organic(f) else "positively_charged_ion"] = 1
            elif q < 0:
                structures_counts["organic_anion" if is_organic(f) else "negatively_charged_ion"] = 1
            else:
                structures_counts["unknown_or_error"] = 1

        # ---- Multiple fragments ----
        else:
            pos_frags = [f for f, q in zip(fragments, frag_charges) if q > 0]
            neg_frags = [f for f, q in zip(fragments, frag_charges) if q < 0]
            neu_frags = [f for f, q in zip(fragments, frag_charges) if q == 0]
            n_solvents = detect_solvents(fragments)

            # Ionic pairs
            if pos_frags and neg_frags and total_charge == 0:
                pos_org = any(is_organic(f) for f in pos_frags)
                neg_org = any(is_organic(f) for f in neg_frags)
                pos_inorg = any(is_transition_metal(f) or not is_organic(f) for f in pos_frags)
                neg_inorg = any(is_transition_metal(f) or not is_organic(f) for f in neg_frags)

                if pos_inorg and neg_inorg:
                    structures_counts["inorganic_salt"] = 1
                elif pos_org and neg_inorg:
                    structures_counts["organic_cation_and_neg_ion"] = 1
                elif neg_org and pos_inorg:
                    structures_counts["organic_anion_and_neg_ion"] = 1
                elif pos_org and neg_org:
                    structures_counts["organic_cation_and_anion"] = 1
                else:
                    structures_counts["anion_and_cation_neutral_pair"] = 1

            elif n_solvents == 1 and len(neu_frags) >= 2:
                structures_counts["molecule_with_solvant"] = 1
            elif n_solvents >= 2 and len(neu_frags) >= 3:
                structures_counts["molecule_with_two_solvants"] = 1
            elif len(neu_frags) >= 3:
                structures_counts["multi_molecule_mixture"] = 1
            else:
                structures_counts["unknown_or_error"] = 1

        # Template
        if '*' in smiles:
            structures_counts["template"] = 1

        return structures_counts

    except Exception:
        return {"unknown_or_error": 1}


@compute_parameter("readiness", "MLIP/NNP/tool readiness", ["uniondict"])
def compute_ToolReadiness(smiles: str,
                          mol: Any,
                          row_raw: Dict[str, object],
                          row_processed: Dict[str, object]):
    """Check whether the molecule can be used with specific MLIPs or NNP tools.
    Returns {tool_name: 1 or 0, optional_reason_keys...} dictionary.
    """

    if mol is None:
        return None

    try:
        readiness: Dict[str, object] = {}

        # Evaluate each tool
        for tool in TOOLS_RULES.keys():
            ok, reason = check_tool_readiness(mol, tool)
            if ok:
                readiness[tool] = 1
        return readiness

    except Exception:
        return {"unknown_or_error": 1}


class ParameterProcessor:
    """
    Process SMILES or RDKit Mol objects through all registered parameter functions.
    """

    def __init__(self, registry: Dict[str, Dict[str, Any]] = PARAMETER_REGISTRY):
        self.registry = registry

    def process(self,
                smiles: str,
                mol: Any,
                row_raw: Dict[str, object],
                row_processed: Dict[str, object],
                trail_dict: Optional[Dict[str, object]] = None) -> Dict[str, object]:
        """
        Calls all registered parameter functions and collects results
        into an ordered dictionary, preserving registration order.
        """
        for name, info in self.registry.items():
            func = info["function"]
            try:
                value = func(smiles, mol, row_raw, row_processed)
            except Exception:
                value = None
            row_processed[name] = value
        return row_processed


def update_aggregates(trail_dict: dict, name: str, value, aggregation: list):
    """Update trail_dict with aggregated values for the given parameter."""
    INIT_SENTINELS = {
        "min": float("inf"),
        "max": float("-inf"),
        "sum": 0,
        "mean": 0,   # mean handled like sum for now; compute final avg elsewhere
        "unionset": set(),
        "uniondict": dict(),
    }

    for agg in aggregation:
        if  agg in {"sum", "mean"}:
            agg = 'sum'
        name_ = f"{name}_{agg}"
        current = trail_dict.get(name_, INIT_SENTINELS[agg])

        if agg == "unionset":
            updated = set(current).union(value if isinstance(value, (list, set)) else [value])
        elif agg == "uniondict":
            updated = dict(current)
            if isinstance(value, (list, set)):# make a shallow copy
                vals = {k:1 for k in value}
            else:
                vals = value

            for k, v in vals.items():
                updated[k] = updated.get(k, 0) + v
        elif agg == "sum":
            updated = current + value
        elif agg == "min":
            updated = min(current, value)
        elif agg == "max":
            updated = max(current, value)
        else:
            raise NotImplementedError(f"Unknown aggregation type: {agg}")

        trail_dict[name_] = updated
    return trail_dict


def sort_elements_universal(
    elements: Union[List[str], Set[str], Dict[str, Any]]
) -> Union[List[str], Dict[str, Any]]:
    """
    Sorts chemical elements (list, set, or dict keys) using RDKit's periodic table order.
    Special symbols ('*', unknowns) are placed first.

    Args:
        elements: list, set, or dict with element symbols (keys if dict)

    Returns:
        Sorted list (for list/set input) or ordered dict (for dict input)
    """
    pt = Chem.GetPeriodicTable()

    def sort_key(el: str):
        # Nonexistent or special symbols → come first
        try:
            atomic_num = pt.GetAtomicNumber(el)
            if el == "*" or atomic_num == 0:
                return (0, el)
            return (1, atomic_num)
        except Exception:
            return (0, el)

    if isinstance(elements, dict):
        sorted_items = sorted(elements.items(), key=lambda kv: sort_key(kv[0]))
        return dict(sorted_items)
    elif isinstance(elements, (set, list)):
        return sorted(list(elements), key=sort_key)
    else:
        raise TypeError("Expected list, set, or dict of element symbols.")

class AggregationProcessor:
    """
    Process SMILES or RDKit Mol objects through all registered parameter functions.
    """
    def __init__(self, registry: Dict[str, Dict[str, Any]] = PARAMETER_REGISTRY):
        self.registry = registry

    def finalize(self, trail_dict: Dict[str, object]):
        for name, info in self.registry.items():
            for agg in info["aggregation"]:
                if agg in {"sum", "mean"}:
                    trail_dict[f"{name}_mean"] = trail_dict[f"{name}_sum"]/trail_dict[f"{name}_len"]

        if 'elements_unionset' in trail_dict:
            trail_dict['elements_unionset'] = sort_elements_universal(trail_dict['elements_unionset'])

        if 'elements_uniondict' in trail_dict:
            trail_dict['elements_uniondict'] = sort_elements_universal(trail_dict['elements_uniondict'])

        if 'elementsinmol_uniondict' in trail_dict:
            trail_dict['elementsinmol_uniondict'] = sort_elements_universal(trail_dict['elementsinmol_uniondict'])

        return trail_dict

    def process(self,
                row_processed: Dict[str, object],
                trail_dict: Dict[str, object]) -> Dict[str, object]:
        """
        Calls all registered parameter functions and collects results
        into an ordered dictionary, preserving registration order.
        """

        for name, info in self.registry.items():
            if row_processed[name] is not None:
                if f'{name}_len' in trail_dict:
                    trail_dict[f'{name}_len'] += 1
                else:
                    trail_dict[f'{name}_len'] = 0
                trail_dict = update_aggregates(trail_dict, name, row_processed[name], info["aggregation"])

        return trail_dict
