from rdkit.Chem import MolFromSmarts, MolFromSmiles, GetFormalCharge
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt


def desalt_building_block(mol):
    def deprotonate_nitrogen(mol):
        """Remove a proton from ammonium species"""
        mol.UpdatePropertyCache()
        patt = MolFromSmarts(
            "[#7+;H1,H2,H3,h1,h2,h3]"
        )  # this pattern matches positive N with at least one proton attached
        try:
            idx = mol.GetSubstructMatches(patt)[0][
                0
            ]  # this raises IndexError if patt is not found
            atom = mol.GetAtomWithIdx(idx)  # get the atom index of the charged N
            atom.SetFormalCharge(0)
            """
            If H are explicit, we have to do explicit removal. If they are implicit, calling UpdatePropertyCache() suffices
            """
            n_hyd = atom.GetNumExplicitHs()
            if n_hyd > 0:
                n_hyd -= 1
                atom.SetNumExplicitHs(n_hyd)
            mol.UpdatePropertyCache()
        except IndexError:
            pass

        return None

    # desalt the building block library
    remover = SaltRemover()
    mol_desalt = remover.StripMol(mol)
    # neutralize ammoniums
    deprotonate_nitrogen(mol_desalt)
    return mol_desalt


def smiles_to_lcms_mass(smiles: str) -> float:
    """
    Calculates the mass of most likely LCMS adduct from the SMILES string of the molecule.

    For neutral molecules, we expect the protonated adduct (M+H)+.
    For positively charged molecules, we expect M+.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        float: The mass of the expected adduct.
    """
    mol = MolFromSmiles(smiles)
    remover = SaltRemover()
    mol_desalt = remover.StripMol(mol)
    charge = GetFormalCharge(mol_desalt)
    if charge == 0:
        return CalcExactMolWt(mol_desalt) + 1.00728  # proton mass
    elif charge == 1:
        return CalcExactMolWt(mol_desalt)
    else:
        raise ValueError(f"Unexpected charge value: {charge}")
