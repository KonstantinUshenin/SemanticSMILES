from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints, rdDetermineBonds
import networkx as nx
import os
import pandas as pd
import argparse
import numpy as np
from rdkit.ML.Cluster import Butina
import tarfile


def mol_to_nx(mol):
    G = nx.Graph()

    # Add nodes (atoms)
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            symbol=atom.GetSymbol(),
            idx=atom.GetIdx(),
        )

    # Add edges (bonds)
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    return G


def check_missing_edges_easy(G_old, G_new):
    """
    Check if G_small is missing edges compared to G_big and if adding these edges would make graphs isomorphic.
    Returns: (is_fixable, missing_edges)
    """
    G_old_edges = set(G_old.edges())
    G_new_edges = set(G_new.edges())
    missing_edges = G_old_edges - G_new_edges
    new_edges = G_new_edges - G_old_edges
    is_new_bonds = bool(new_edges)
    is_broken_bonds = bool(missing_edges)
    return is_new_bonds, is_broken_bonds, list(new_edges), list(missing_edges)


def gen_conformers(
    mol,
    numConfs=100,
    maxAttempts=10,
    pruneRmsThresh=0.1,
    useExpTorsionAnglePrefs=True,
    useBasicKnowledge=True,
    enforceChirality=True,
):
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xF00D
    ps.numThreads = 4
    ps.maxAttempts = maxAttempts
    ps.pruneRmsThresh = pruneRmsThresh
    ps.useExpTorsionAnglePrefs = useExpTorsionAnglePrefs
    ps.useBasicKnowledge = useBasicKnowledge
    ps.enforceChirality = enforceChirality
    ids = AllChem.EmbedMultipleConfs(mol, numConfs, ps)
    return list(ids)


def write_conformers_to_sdf(mol, filename, rmsClusters, conformerPropsDict, minEnergy):
    w = Chem.SDWriter(filename)
    for cluster in rmsClusters:
        for confId in cluster:
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            conformerProps = conformerPropsDict[confId]
            mol.SetIntProp("conformer_id", confId + 1)
            for key in conformerProps.keys():
                mol.SetProp(key, str(conformerProps[key]))
            e = conformerProps["energy_abs"]
            if e:
                mol.SetDoubleProp("energy_delta", e - minEnergy)
            w.write(mol, confId=confId)
    w.flush()
    w.close()


def calc_energy(mol, conformerId, minimizeIts):
    ff = AllChem.MMFFGetMoleculeForceField(
        mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId
    )
    ff.Initialize()
    ff.CalcEnergy()
    results = {}
    if minimizeIts > 0:
        results["converged"] = ff.Minimize(maxIts=minimizeIts)
    results["energy_abs"] = ff.CalcEnergy()
    return results


def cluster_conformers(mol, mode="RMSD", threshold=2.0):
    if mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(mol)
    else:
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    rms_clusters = Butina.ClusterData(
        dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True
    )
    return rms_clusters


def align_conformers(mol, clust_ids):
    rmslist = []
    AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
    return rmslist


def add_file_to_archive(filename, archivename):
    if os.path.exists(archivename):
        mode = "a"
    else:
        mode = "w"
    with tarfile.open(archivename, mode=mode) as tf:
        tf.add(filename)
    os.remove(filename)


def write_conformers_to_xyz(mol, filename, rmsClusters):
    with open(filename, "a") as w:
        for cluster in rmsClusters:
            output = Chem.rdmolfiles.MolToXYZBlock(mol, confId=cluster[0])
            w.write(output)


def generate_conformers(ds, split, archivename="centroids.tar"):
    numConfs = 100
    maxAttempts = 10
    pruneRmsThresh = 0.1
    clusterMethod = "RMSD"
    clusterThreshold = 2.0
    minimizeIterations = 0
    last_computed = -1
    if os.path.exists(archivename):
        last_computed = tarfile.open(archivename).getmembers()[-1].name
        last_computed = int(last_computed.split("_")[0])
    # print (last_computed)
    i = 0
    for sm_id in tqdm(split):
        if sm_id <= last_computed + 1:
            continue

        smile = ds["SMILES"][sm_id]
        # print (smile)
        i = i + 1
        mol = Chem.MolFromSmiles(smile)
        m = Chem.AddHs(mol)
        # generate the confomers
        conformerIds = gen_conformers(
            m, numConfs, maxAttempts, pruneRmsThresh, True, True, True
        )
        if not conformerIds:
            continue

        removed = 0
        currentId = 0
        for confId in range(m.GetNumConformers()):
            mol_xyz = Chem.Mol()
            editable_mol = Chem.EditableMol(mol_xyz)

            # Add atoms with their 3D positions
            for atom in m.GetAtoms():
                atom_idx = editable_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))

            mol_xyz = editable_mol.GetMol()
            conf = Chem.Conformer(m.GetNumAtoms())
            for i, pos in enumerate(m.GetConformer(confId).GetPositions()):
                conf.SetAtomPosition(i, pos)
            mol_xyz.AddConformer(conf)
            rdDetermineBonds.DetermineConnectivity(mol_xyz)
            mol_smi_nx = mol_to_nx(m)
            mol_xyz_nx = mol_to_nx(mol_xyz)
            is_new_bonds, is_broken_bonds, new_edges, missing_edges = (
                check_missing_edges_easy(mol_smi_nx, mol_xyz_nx)
            )
            if is_new_bonds or is_broken_bonds:
                m.RemoveConformer(confId)
                removed += 1
            else:
                m.GetConformer(confId).SetId(currentId)
                currentId += 1
        # print (removed)
        if m.GetNumConformers() == 0:
            continue

        rmsClusters = cluster_conformers(m, clusterMethod, clusterThreshold)

        clustersSizes = [len(cluster) for cluster in rmsClusters]

        # print ("Molecule", i, ": generated", len(conformerIds), "conformers and", len(rmsClusters), "clusters")

        clustersSizes = [len(cluster) for cluster in rmsClusters]
        max_id = np.argmax(np.cumsum(clustersSizes) / m.GetNumConformers() > 0.95)
        rmsClustersPerCluster = []
        clusterNumber = 0
        for cluster in rmsClusters:
            clusterNumber = clusterNumber + 1
            rmsWithinCluster = align_conformers(m, cluster)

        write_conformers_to_xyz(
            m, str(sm_id) + "_centroids.xyz", rmsClusters[: max_id + 1]
        )
        add_file_to_archive(str(sm_id) + "_centroids.xyz", archivename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make conformers")
    parser.add_argument("b", type=int, help="batch_number")
    args = parser.parse_args()
    ds = pd.read_csv("/mnt/sdd/khrabrov/zinc/zinc_moseslike_newats_250K.csv")

    generate_conformers(
        ds,
        range(args.b * 50000, (args.b + 1) * 50000),
        f"/mnt/sdd/data/zinc_conformers/archive_{args.b}.tar",
    )
