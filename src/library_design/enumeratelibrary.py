"""
Enumerate combinatorial libraries and dump products to SDF.
(This substitutes the Schrödinger workflow)

Inputs:
    - library_constituents_dataframe.pkl: DataFrame of rdkit mol objects, output of inventorytosdf.py

Outputs:
    - products_[A-H].sdf: SDF of products for all reagent combinations, sorted by type where A is the main (expected)
      product and all others are side products

Actions:
    - Read pickle of DataFrame
    - Desalt + neutralize molecules. Neutralization only affects nitrogen cations
    - Define reactions from reactionSMARTS
    - Check if all input molecules are valid reactants + no duplicates
    - Enumerate library
    - Sanitize + rearrange products for output
    - Output to SDF
"""

import pandas as pd
import gzip
import warnings
import string
import itertools
import pickle as pkl

from rdkit import Chem
from rdkit.Chem import Draw, SaltRemover, AllChem
from rdkit.Chem.SimpleEnum.Enumerator import EnumerateReaction

from src.definitions import LIB_INFO_DIR, LIB_SDF_DIR

"""GLOBALS"""
POSTPROCESSING_ONLY = True
DEBUG = 1


def deprotonate_nitrogen(mol):
    """Remove a proton from ammonium species"""
    mol.UpdatePropertyCache()
    if DEBUG >= 2:
        print(f"Next: {AllChem.CalcMolFormula(mol)}")
    patt = Chem.MolFromSmarts(
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
        if DEBUG >= 2:
            print("No charged nitrogen found")
    if DEBUG >= 2:
        print(f".....{AllChem.CalcMolFormula(mol)}")
    return None


def check_reactants(rxn, name, *reactants):
    """
    Check if any of the desired reactants does not work as a reactant in a reaction
    :param rxn: rdkit reaction object
    :param name: str
    :param reactants: list of pd.DataFrames
    :return passed: bool
    """
    r = pd.concat(reactants)
    not_reactant = r[~r.loc[:, "desalted"].apply(rxn.IsMoleculeReactant)]
    print(r["Compound Name"])
    if len(not_reactant) > 0:
        passed = False
        warnings.warn(
            "A molecule could not be assigned as reactant with the current reaction definition.",
            Warning,
        )
        print(f"Reaction: {name}")
        print("Unassigned reactants:")
        print(not_reactant["Compound Name"])
    else:
        passed = True
        print(f"Reactant check passed for reaction {name}")
    return passed


def check_unique_molecules(mols: list):
    uniq_mol = {}
    duplicate_mol = []
    for i in mols:
        smi = Chem.MolToSmiles(i)
        if smi in uniq_mol.keys():
            warnings.warn("A building block has a duplicate:", Warning)
            print(smi)
            duplicate_mol.append(i)
        uniq_mol[smi] = i
    return duplicate_mol


def add_name_prop_to_mol(reactant1, reactant2, reactant3, product_library):
    """
    Add name properties to enumerated molecules. The cartesian enumeration strategy is (0,0,0), (1,0,0),... according
    to the documentation. In reality I observe that the strategy is (0,0,0), (0,0,1),... which is used for assigning
    names in this function
    """
    names = [
        f'{i.GetProp("_Name")} + {m.GetProp("_Name")} + {t.GetProp("_Name")}'
        for i in reactant1
        for m in reactant2
        for t in reactant3
    ]  # generate the names from reactant properties
    if len(product_library) != len(
        names
    ):  # check if the correct number of names was generated
        raise ValueError(
            f"Number of products ({len(product_library)}) does not fit number of inferred names ({len(names)})"
        )
    [
        p.SetProp("_Name", n)
        for p_list, n in zip(product_library, names)
        for p in p_list
    ]  # assign the name property
    return product_library


if POSTPROCESSING_ONLY is False:
    """Import from pickle"""
    compounds = pd.read_pickle(LIB_INFO_DIR / "library_constituents_dataframe.pkl")

    """Desalt building blocks and deprotonate N"""
    # desalt the building block library
    remover = SaltRemover.SaltRemover()
    compounds["desalted"] = compounds.loc[:, "mol"].apply(remover.StripMol)
    # neutralize ammoniums
    compounds.loc[:, "desalted"].apply(deprotonate_nitrogen)
    compounds["desalted_SMILES"] = compounds.loc[:, "desalted"].apply(Chem.MolToSmiles)

    """Define reactions"""
    # all 8 TerTH prods
    rxn_TH = AllChem.ReactionFromSmarts(
        "F[B-](F)(F)[#6](-[*:1])=O.O=[#6]1-[#8]C2([#6]-[#6]-[#6]-[#6]-[#6]2)[#8]C11[#6:3]-[#6:2]-[#7]-[#8]1.[#6:4]-[#6](=S)-[#7]-[#7]>>[#6:4]-[#6]-1=[#7]-[#7]=[#6](-[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O)-[#16]-1.[#6:4]-[#6]-1=[#7]-[#7]C([#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O)([#16]-1)[#6](-[#8])=O.[#6:4]-[#6]-1=[#7]-[#7+]2=[#6](-[*:1])-[#7]-[#6:2]-[#6:3]C2([#16]-1)[#6](-[#8-])=O.[#6:4]-[#6]-1=[#7]-[#7]=[#6](-[*:1])-[#16]-1.[#6:4]-[#6]-1=[#7]-[#7]=[#6](-[#6:4])-[#16]-1.[#8]-[#6](=O)-[#6](=O)-[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O.[#8]-[#6](=O)-[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O.[#6:4]-[#6]-1=[#7]-[#7]=[#6](-[#16]-1)-[#6:3]=[#6:2]"
    )

    # all 8 TerABT prods
    rxn_ABT = AllChem.ReactionFromSmarts(
        "F[B-](F)(F)[#6](-[*:1])=O.O=[#6]1-[#8]C2([#6]-[#6]-[#6]-[#6]-[#6]2)[#8]C11[#6:3]-[#6:2]-[#7]-[#8]1.[#7]-c1[c:4][c:5][c:6][c:7]c1-[#16]>>[*:1]-[#6](=O)-[#7]-[#6:2]-[#6:3]-c1nc2[c:4][c:5][c:6][c:7]c2s1.[#8]-[#6](=O)C1([#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O)[#7]-c2[c:4][c:5][c:6][c:7]c2-[#16]1.[#8-]-[#6](=O)C12[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=[#7+]1-c1[c:4][c:5][c:6][c:7]c1-[#16]2.[*:1]-c1nc2[c:4][c:5][c:6][c:7]c2s1.[#7]-c1[c:4][c:5][c:6][c:7]c1-[#16]-[#16]-c1[c:7][c:6][c:5][c:4]c1-[#7].[#8]-[#6](=O)-[#6](=O)-[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O.[#8]-[#6](=O)-[#6:3]-[#6:2]-[#7]-[#6](-[*:1])=O.[#6:2]=[#6;v4:3]-c1nc2[c:4][c:5][c:6][c:7]c2s1"
    )

    # prepare for visualization
    AllChem.Compute2DCoordsForReaction(rxn_TH)
    AllChem.Compute2DCoordsForReaction(rxn_ABT)
    # prepare for enumeration
    rxn_TH.Initialize()
    rxn_ABT.Initialize()
    n_warn_TH, n_err_TH = rxn_TH.Validate(silent=True)
    n_warn_ABT, n_err_ABT = rxn_ABT.Validate(silent=True)
    if n_err_TH > 0:
        raise ValueError(f"Invalid reaction gave {n_err_TH} errors in validation")
    if n_err_ABT > 0:
        raise ValueError(f"Invalid reaction gave {n_err_ABT} errors in validation")

    """Control reactions visually"""
    Draw.ReactionToImage(rxn_TH)
    Draw.ReactionToImage(rxn_ABT)

    """define reagents"""
    all_KAT = compounds[compounds.loc[:, "Category"].str.startswith("I")]
    all_Mon = compounds[compounds.loc[:, "Category"].str.startswith("M")]
    all_Spiro = compounds[compounds.loc[:, "Compound Name"].str.startswith("Spiro")]
    all_Fused = compounds[compounds.loc[:, "Compound Name"].str.startswith("Fused")]
    all_Sub = compounds[compounds.loc[:, "Compound Name"].str.startswith("Mon")]
    all_TerTH = compounds[compounds.loc[:, "Compound Name"].str.startswith("TerTH")]
    all_TerABT = compounds[compounds.loc[:, "Compound Name"].str.startswith("TerABT")]

    I = all_KAT["desalted"].tolist()
    M = all_Mon["desalted"].tolist()
    fused = all_Fused["desalted"].tolist()
    spiro = all_Spiro["desalted"].tolist()
    sub = all_Sub["desalted"].tolist()
    T_TH = all_TerTH["desalted"].tolist()
    T_ABT = all_TerABT["desalted"].tolist()

    """print reagent info"""
    print(f"KATs: {len(I)}")
    print(f"Mons: {len(M)}")
    print(f" - Spiro: {len(spiro)}")
    print(f" - Fused: {len(fused)}")
    print(f" - Sub: {len(sub)}")
    print(f"TerTHs: {len(T_TH)}")
    print(f"TerABTs: {len(T_ABT)}\n")
    print(f"expected TH (A) products: {len(I) * len(M) * len(T_TH)}")
    print(f" - expected TH (A) products (spiro): {len(I) * len(spiro) * len(T_TH)}")
    print(f" - expected TH (A) products (fused): {len(I) * len(fused) * len(T_TH)}")
    print(f" - expected TH (A) products (sub): {len(I) * len(sub) * len(T_TH)}\n")
    print(f"expected ABT (A) products: {len(I) * len(M) * len(T_ABT)}")
    print(f" - expected ABT (A) products (spiro): {len(I) * len(spiro) * len(T_ABT)}")
    print(f" - expected ABT (A) products (fused): {len(I) * len(fused) * len(T_ABT)}")
    print(f" - expected ABT (A) products (sub): {len(I) * len(sub) * len(T_ABT)}")

    check_reactants(rxn_TH, "TH", all_KAT, all_Mon, all_TerTH)
    check_reactants(rxn_ABT, "ABT", all_KAT, all_Mon, all_TerABT)

    """check uniqueness of building blocks"""
    duplicate_I = check_unique_molecules(I)
    duplicate_M = check_unique_molecules(M)
    duplicate_T = check_unique_molecules(T_TH + T_ABT)

    """run enumeration (consumes lots of RAM)"""
    reactant_I = I
    reactant_M = M
    reactant_ABT = T_ABT
    reactant_TH = T_TH
    product_generator_ABT = EnumerateReaction(
        rxn_ABT, (reactant_I, reactant_M, reactant_ABT), uniqueProductsOnly=True
    )
    product_generator_TH = EnumerateReaction(
        rxn_TH, (reactant_I, reactant_M, reactant_TH), uniqueProductsOnly=True
    )
    product_list_ABT = list(product_generator_ABT)
    product_list_TH = list(product_generator_TH)

    product_list_ABT = add_name_prop_to_mol(
        reactant_I, reactant_M, reactant_ABT, product_list_ABT
    )
    product_list_TH = add_name_prop_to_mol(
        reactant_I, reactant_M, reactant_TH, product_list_TH
    )
    product_list = product_list_ABT + product_list_TH

    """Pickle the library, just in case postprocessing fails."""
    with open(LIB_SDF_DIR / "product_list_unprocessed.pkl", "wb") as file:
        pkl.dump(product_list, file)

"""Reuse the pickled library if the option is set"""
if POSTPROCESSING_ONLY is True:
    print("Loading product list from pickle...")
    with open(LIB_SDF_DIR / "product_list_unprocessed.pkl", "rb") as file:
        product_list = pkl.load(file)
    print("Finished loading.")

"""Post-process the generated library"""
# sanitize mols. We don't use a list comprehension to be able to catch exceptions
errors = []
for i, p_list in enumerate(product_list):
    for j, p in enumerate(p_list):
        problems = Chem.DetectChemistryProblems(p)
        if len(problems) > 0:
            if problems[0].GetType() == "AtomValenceException":
                atom_idx = problems[0].GetAtomIdx()
                p.GetAtomWithIdx(atom_idx).SetNumExplicitHs(0)
        try:
            Chem.SanitizeMol(p)
        except:
            # if Sanitization throws exception, replace the molecule with water.
            # we need to take a little detour as p_list is a tuple.
            temp_list = list(p_list)
            temp_list[j] = Chem.MolFromSmiles("O")  # water acts as a placeholder
            product_list[i] = tuple(temp_list)
            # product_list[i] = (p_list[0],
            #                    p_list[1],
            #                    p_list[2],
            #                    p_list[3],
            #                    p_list[4],
            #                    p_list[5],
            #                    p_list[6],
            #                    Chem.MolFromSmiles('O'),  # this is the one we want to replace
            #                    )

resorted_products = list(
    map(list, itertools.zip_longest(*product_list, fillvalue=None))
)  # "transposes" the list of lists


"""Write sdf files. One file per product type."""
for p, s in zip(resorted_products, string.ascii_uppercase):
    print(f"Saving products {s}...")
    with gzip.open(LIB_SDF_DIR / f"product_{s}.sdf.gz", "wt") as file:
        writer = Chem.SDWriter(file)
        for mol in p:
            writer.write(mol)
        writer.close()

    print(f"Products {s} were saved to SDF")
