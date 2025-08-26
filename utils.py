from Bio import PDB
import numpy as np
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Residue import Residue
from collections import defaultdict

# parser = PDB.PDBParser(QUIET=True)
parser = PDB.MMCIFParser(QUIET=True)
# Atomic radii for various atom types.
# You can comment out the ones you don't care about or add new ones
atom_radii = {
    #    "H": 1.20,  # Who cares about hydrogen??
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "F": 1.47,
    "P": 1.80,
    "CL": 1.75,
    "MG": 1.73,
}
_WATER = {"HOH", "WAT", "DOD", "H2O"}


def _is_protein_residue(res: Residue) -> bool:
    """Standard amino-acid residue: hetflag == ' ' and is_aa True"""
    return res.id[0] == " " and is_aa(res, standard=True)


def _is_ligand_residue(
    res: Residue,
    ligand_resnames: Collection[str] | None = None,
    ignore_waters: bool = True,
) -> bool:
    """Check if residue is from a ligand"""
    resname = res.resname.strip()
    hetflag = res.id[0]
    if ligand_resnames is not None:
        return resname in ligand_resnames
    # Default: any hetero (HETATM) that isn't water
    if hetflag != " ":
        if ignore_waters and resname in _WATER:
            return False
        return True
    return False


def unique_residues_per_ligand(clashes_dict):
    per_lig = defaultdict(set)
    for (lig_id, _), hits in clashes_dict.items():
        per_lig[lig_id].update(resid[1] for _, resid, _ in hits)
    return {lig: sorted(v) for lig, v in per_lig.items()}


def count_clashes(structure, clash_cutoff=0.63):
    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {
        i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j]))
        for i in atom_radii
        for j in atom_radii
    }
    # Extract protein and ligand atoms for which we have a radii
    protein_atoms = [
        x
        for x in structure.get_atoms()
        if x.element in atom_radii and _is_protein_residue(x.get_parent())
    ]
    protein_coords = np.array([a.coord for a in protein_atoms], dtype="d")
    ligand_atoms = [
        x
        for x in structure.get_atoms()
        if x.element in atom_radii and _is_ligand_residue(x.get_parent())
    ]
    ligand_coords = np.array([a.coord for a in ligand_atoms], dtype="d")

    # Build a KDTree for protein atoms (speedy!!!)
    kdt = PDB.kdtrees.KDTree(protein_coords)

    # Initialize a dict to hold clashes
    clashes_dict = {}
    # Iterate through all atoms
    for lig_atm in ligand_atoms:
        thresholds = [
            (k, v) for k, v in clash_cutoffs.items() if k[0] == lig_atm.element
        ]
        for pair, t in thresholds:
            # Find atoms that could be clashing
            kdt_search = kdt.search(np.array(lig_atm.coord, dtype="d"), t)
            clashes = [
                (
                    protein_atoms[a.index].get_serial_number(),
                    protein_atoms[a.index].get_parent().get_id(),
                    protein_atoms[a.index].element,
                )
                for a in kdt_search
                if protein_atoms[a.index].element == pair[-1]
            ]
            if len(clashes) > 0:
                clashes_dict[(lig_atm.get_parent().get_id()[0], lig_atm.get_id())] = (
                    clashes
                )

    return unique_residues_per_ligand(clashes_dict)


structure = parser.get_structure("model", "../../../Downloads/clash_test.cif")
# structure = parser.get_structure("model", "../../../Downloads/MAPK_poses_maxit/0i54/0i54.pdb")
clashes = count_clashes(structure)
print(clashes)
