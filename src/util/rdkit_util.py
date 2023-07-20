from rdkit import Chem
from rdkit.Chem import MolFromSmarts, MolFromSmiles, GetFormalCharge, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
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
    elif charge == -1:
        return CalcExactMolWt(mol_desalt) + 2 * 1.00728
    else:
        raise ValueError(f"Unexpected charge value: {charge}")


def create_reaction_instance(rxn, reactants):
    """
    Create an instance of a reaction, given reactants, and map all atoms that end up in the product(s).
    This is adapted from Greg's code in https://github.com/rdkit/rdkit/issues/1269#issuecomment-275228746,
    but extended to map the side chains as well.
    Copied from https://github.com/jugoetz/slap-platform-predict/blob/5acb77ff4fdc412f7fd03a8226a4635e086978d7/src/util/rdkit_util.py#L82
    Note that atoms that are not present in the product (unbalanced reaction equation) will not be annotated.
    """

    # first, we set a tag on reactant atoms. This will be passed on to the product for all non-mapped atoms
    for i, sm in enumerate(reactants):
        for atom in sm.GetAtoms():
            atom.SetProp("tag", "reactant-%s atom-%s" % (i, atom.GetIdx()))

    # for the mapped atoms, extract their mapping in the reactants
    map_number_to_reactant = {}
    for i, reactant in enumerate(rxn.GetReactants()):
        for atom in reactant.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_number_to_reactant[atom.GetIntProp("molAtomMapNumber")] = (
                    i,
                    atom.GetIdx(),
                )

    mapped_reactions = []  # this will hold the reactions
    product_set = rxn.RunReactants(reactants)  # run the reaction to get product set
    if len(product_set) == 2:
        # if one of the reagents is cyclic?or a ketone? we get two identical product sets. We remove one of them.
        if Chem.MolToSmiles(product_set[0][0]) == Chem.MolToSmiles(product_set[1][0]):
            product_set = [
                product_set[0],
            ]

    # now, we look into the products
    for products in product_set:
        # we need to know the highest mapno, because mapping the "tagged" atoms will have to start above that
        mapno_max = max(
            map_number_to_reactant.keys()
        )  # needs to reset on every product_set
        reactant_list = [Chem.Mol(x) for x in reactants]
        reaction = rdChemReactions.ChemicalReaction()
        for p in products:
            for atom in p.GetAtoms():
                # for atoms that are mapped in the reaction template
                # (RDKit does not copy user-defined properties to the products for mapped atoms, only for unmapped atoms)
                # the reference solution from Greg uses "if not atom.HasProp('old_mapno'):" here, but that is too wide
                # and fails for cyclic ketones, where the carbonyl C will be mapped twice.
                if not atom.HasProp("tag"):
                    mno = atom.GetIntProp("old_mapno")
                    atom.SetIntProp("molAtomMapNumber", mno)
                    ridx, aidx = map_number_to_reactant[mno]
                    # aidx is the index of the atom in the reactant template. We need
                    # to read out the number in the actual reactant:
                    raidx = atom.GetIntProp("react_atom_idx")
                    ratom = (
                        reactant_list[ridx]
                        .GetAtomWithIdx(raidx)
                        .SetIntProp("molAtomMapNumber", mno)
                    )

                # for atoms that are unmapped in the reaction template
                else:
                    tag = atom.GetProp("tag")
                    mapno_max += 1
                    atom.SetIntProp("molAtomMapNumber", mapno_max)
                    # now find the tag in reactant_list
                    for sm in reactant_list:
                        for ratom in sm.GetAtoms():
                            if ratom.HasProp("tag"):
                                if ratom.GetProp("tag") == tag:
                                    ratom.SetIntProp("molAtomMapNumber", mapno_max)

            # now add the product to the reaction
            reaction.AddProductTemplate(p)
        # add the reactants to reaction
        for reactant in reactant_list:
            reaction.AddReactantTemplate(reactant)
        # add reaction for all product sets
        mapped_reactions.append(reaction)
    return mapped_reactions


def map_reactions(rxn, reactant_sets):
    """
    Take a reaction template and a list of reactant sets and return the mapped reactions.
    Adapted from https://github.com/jugoetz/slap-platform-predict/blob/5acb77ff4fdc412f7fd03a8226a4635e086978d7/src/util/rdkit_util.py#L163
    """
    mapped_reactions = []
    for i, reactant_set in enumerate(reactant_sets):
        reaction_inst = create_reaction_instance(rxn, reactant_set)
        if len(reaction_inst) == 1:  # all good
            mapped_reactions.append(reaction_inst[0])
        elif len(reaction_inst) == 0:  # failed
            mapped_reactions.append(None)
            print(f"ERROR: No product for reactant set with index {i}")
        else:  # too many resulting reactions
            # remove any duplicates
            idx = []
            unique_reactions = []
            for j, reac in enumerate(reaction_inst):
                reac_smarts = rdChemReactions.ReactionToSmarts(reac)
                if reac_smarts not in unique_reactions:
                    idx.append(j)
                    unique_reactions.append(reac_smarts)
                reaction_inst_cleaned = [reaction_inst[j] for j in idx]
                if len(reaction_inst_cleaned) > 1:
                    print(f"ERROR: Multiple products for reactant set with index {i}")

            mapped_reactions.append(reaction_inst_cleaned)

    return mapped_reactions


def draw_chemical_reaction(smiles, highlightByReactant=True, font_scale=1.5):
    """
    Draw reactions in a nice fashion with atom map numbers and optional reactant highlighting

    https://gist.github.com/greglandrum/61c1e751b453c623838759609dc41ef1
    """

    def moveAtomMapsToNotes(m):
        """Move atom maps to be annotations (so they can be drawn)"""
        for at in m.GetAtoms():
            if at.GetAtomMapNum():
                at.SetProp("atomNote", str(at.GetAtomMapNum()))

    rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
    trxn = rdChemReactions.ChemicalReaction(rxn)
    for m in trxn.GetReactants():
        moveAtomMapsToNotes(m)
    for m in trxn.GetProducts():
        moveAtomMapsToNotes(m)
    d2d = rdMolDraw2D.MolDraw2DSVG(800, 300)
    d2d.drawOptions().annotationFontScale = font_scale
    d2d.DrawReaction(trxn, highlightByReactant=highlightByReactant)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()
