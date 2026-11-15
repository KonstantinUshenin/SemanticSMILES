from typing import List
import pandas as pd
import re

import json

try:
    # Package-style import (e.g., `import notebooks.processing`)
    from .analysis import ParameterProcessor, AggregationProcessor
except ImportError:
    # Notebook-folder import (e.g., cwd=notebooks, `import processing`)
    from analysis import ParameterProcessor, AggregationProcessor


def split_smiles(smiles: str):
    """
    Split a possibly multi-molecule or reaction SMILES string into
    a list of canonical molecule SMILES.

    Handles:
      - Dot-separated mixtures/salts/solvates
      - Reaction SMILES (>, >>, ->)
      - Nonstandard delimiters (;, |, ,)
      - Ionic fragments (e.g., [Na+].[Cl-])
      - Cleans and canonicalizes unique molecules
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return []

    # --- Normalize separators ---
    s = smiles.strip()
    # Replace reaction arrows
    s = re.sub(r">+|->", ".", s)
    # Replace nonstandard delimiters
    s = re.sub(r"[;|,]", ".", s)

    # --- Split by dot into fragments ---
    fragments = [frag.strip() for frag in s.split(".") if frag.strip()]

    mols = []
    for frag in fragments:
        try:
            mol = Chem.MolFromSmiles(frag)
            if mol:
                # Canonicalize
                mols.append(Chem.MolToSmiles(mol))
        except Exception:
            continue  # skip invalid

    # --- Deduplicate and sort for stability ---
    return sorted(set(mols))


from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import inchi
import logging
from io import StringIO

def normalize_smiles_and_inchikey(smiles: str):
    """
    Normalize SMILES and compute InChIKey, capturing RDKit log messages.

    Returns:
        (norm_smiles | None,
         inchikey | None,
         has_warning: bool,
         log_text: str | None,
         mol_with_h: Chem.Mol | None)
    """
    if not smiles or not isinstance(smiles, str):
        return None, None, True, "Invalid input type or empty string", None

    rdBase.LogToPythonLogger()
    rdkit_logger = logging.getLogger("rdkit")

    stream = StringIO()
    handler = logging.StreamHandler(stream)
    handler.setLevel(logging.INFO)
    fmt = logging.Formatter("%(levelname)s: %(message)s")
    handler.setFormatter(fmt)

    old_level = rdkit_logger.level
    rdkit_logger.setLevel(logging.INFO)
    rdkit_logger.addHandler(handler)

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log_text = stream.getvalue().strip() or "RDKit could not parse SMILES"
            return None, None, True, log_text, None

        # Add hydrogens
        norm_smiles = Chem.MolToSmiles(mol, canonical=True)
        inchikey = inchi.MolToInchiKey(mol)
        mol_h = Chem.AddHs(mol)

        log_text = stream.getvalue().strip() or None
        has_warning = bool(log_text)
        return norm_smiles, inchikey, has_warning, log_text, mol_h

    except Exception as e:
        log_text = stream.getvalue().strip()
        msg = f"{e}" if not log_text else f"{log_text}\n{e}"
        return None, None, True, msg, None

    finally:
        rdkit_logger.removeHandler(handler)
        rdkit_logger.setLevel(old_level)
        handler.flush()
        handler.close()


def iter_smiles_from_chunk(chunk, smiles_columns_list):
    """
    Stream SMILES strings from a large CSV file in chunks.
    Yields tuples: (row_index, column_name, smiles_string)
    """
    for idx, row in chunk.iterrows():
        for col in smiles_columns_list:
            raw_smiles = row[col]
            smiles_list = split_smiles(raw_smiles)  # your existing function
            for smiles in smiles_list:
                yield idx, col, smiles


class BasicAnalysis:
    """Minimal analysis of data"""

    def __init__(self,
                 dataset_key,
                 file_in:str,
                 smiles_columns_list: List[str],
                 file_full_desc:str,
                 file_unique_desc:str,
                 file_broken_desc:str,
                 file_dataset_summary:str,
                 chunk_size=50000) -> None:
        self.dataset_key = dataset_key
        self.file_in = file_in
        self.smiles_columns_list = smiles_columns_list
        self.file_full_desc = file_full_desc
        self.file_unique_desc = file_unique_desc
        self.file_broken_desc = file_broken_desc
        self.file_dataset_summary = file_dataset_summary
        self.chunk_size = chunk_size

        self.processor = ParameterProcessor()
        self.aggregation = AggregationProcessor()

    def process(self):

        GLOBAL_DICT_OF_UNIQUE = {None:{None:0}} # InChIKey: {SMILES:None}

        GLOBAL_STATISTICS_FULL = {"dataset_key":self.dataset_key,
                                  "number_of_records":0,
                                  "number_of_molecules": 0,
                                  "number_of_unique_molecules": 0,
                                  "number_of_SMILES_with_error":0,
                                  "number_of_SMILES_with_warning":0}

        # GLOBAL_STATISTICS_UNIQUE = {"number_of_rows":0,
        #                             "number_of_molecules": 0,
        #                             "number_of_broken_SMILES":0,
        #                             "number_of_incorrect_SMILES":0}

        _header_written_broken_SMILES = False
        _header_written_full_SMILES = False
        _header_written_unique_SMILES = False

        # Processing each unique molecules from all columns and all SMILES in a separate mode
        for chunk in pd.read_csv(self.file_in, chunksize=self.chunk_size):
            GLOBAL_STATISTICS_FULL['number_of_records'] += len(chunk)

            list_dict_broken = []
            list_dict_full = []
            list_dict_unique = []
            for idx, col_orig, smiles in iter_smiles_from_chunk(chunk, self.smiles_columns_list):
                norm_SMILES, InChIKey, is_warn, msg, mol_h = normalize_smiles_and_inchikey(smiles)

                # Store broken smiles to dict
                if (norm_SMILES is None) or is_warn:
                    dict_broken = {"SMILES": smiles,
                                   "norm_SMILES": norm_SMILES,
                                   "is_broken": (norm_SMILES is None),
                                   "is_warning": is_warn,
                                   "msg": msg}
                    list_dict_broken.append(dict_broken)
                    if (norm_SMILES is None):
                        GLOBAL_STATISTICS_FULL["number_of_SMILES_with_error"] += 1
                    else:
                        GLOBAL_STATISTICS_FULL["number_of_SMILES_with_warning"] += 1


                if norm_SMILES is not None:
                    # Store NON-UNIQUE molecules to dict
                    dict_full = {"norm_SMILES": norm_SMILES,
                                 "InChIKey": InChIKey}
                    dict_full = self.processor.process(norm_SMILES, mol_h, col_orig, dict_full)
                    list_dict_full.append(dict_full)
                    GLOBAL_STATISTICS_FULL["number_of_molecules"] += 1

                    # Store UNIQUE SMILES to dict
                    _is_store = True
                    if InChIKey in GLOBAL_DICT_OF_UNIQUE:
                        _is_store = not (norm_SMILES in GLOBAL_DICT_OF_UNIQUE[InChIKey])
                    else:
                        _is_store = True

                    if _is_store:
                        list_dict_unique.append(dict_full)
                        GLOBAL_STATISTICS_FULL["number_of_unique_molecules"] += 1
                        GLOBAL_STATISTICS_FULL = self.aggregation.process(dict_full, GLOBAL_STATISTICS_FULL)

                        if InChIKey not in GLOBAL_DICT_OF_UNIQUE:
                            GLOBAL_DICT_OF_UNIQUE[InChIKey] = {norm_SMILES:0}
                        else:
                            if norm_SMILES not in GLOBAL_DICT_OF_UNIQUE[InChIKey]:
                                GLOBAL_DICT_OF_UNIQUE[InChIKey][norm_SMILES] = 0
                            else:
                                GLOBAL_DICT_OF_UNIQUE[InChIKey][norm_SMILES] += 1

            # Store broken smiles to file
            if len(list_dict_broken) > 0:
                df_processed = pd.DataFrame(list_dict_broken)
                df_processed['dataset_key'] = self.dataset_key
                df_processed.to_csv(
                    self.file_broken_desc,
                    mode="a",
                    index=False,
                    header=not _header_written_broken_SMILES,
                )
                _header_written_broken_SMILES = True

            # Store unique smiles to file
            if len(list_dict_unique) > 0:
                df_processed = pd.DataFrame(list_dict_unique)
                df_processed['dataset_key'] = self.dataset_key
                df_processed.to_csv(
                    self.file_unique_desc,
                    mode="a",
                    index=False,
                    header=not _header_written_unique_SMILES,
                )
                _header_written_unique_SMILES = True

            # Store full smiles to file
            if len(list_dict_full)>0:
                df_processed = pd.DataFrame(list_dict_full)
                df_processed['dataset_key'] = self.dataset_key
                df_processed.to_csv(
                    self.file_full_desc,
                    mode="a",
                    index=False,
                    header=not _header_written_full_SMILES,
                )
                _header_written_full_SMILES = True

            # Free memory
            del list_dict_broken
            del list_dict_full
            del list_dict_unique

        GLOBAL_STATISTICS_FULL = self.aggregation.finalize(GLOBAL_STATISTICS_FULL)
        with open(self.file_dataset_summary, "w") as json_file:
            json.dump(
                GLOBAL_STATISTICS_FULL,
                json_file,
                default=lambda o: list(o) if isinstance(o, set) else o,
                indent=2
            )
