import copy
import re
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import MolFromSmarts


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


def parse_formula(formula: str) -> dict:  # Formula Parsing by Aditya Matam
    def multiply(formula: dict, mul: int) -> None:
        for key in formula:
            formula[key] *= mul

    formDict = {}
    # PARENS
    for match in re.finditer(r"\((.*?)\)(\d*)", formula):
        parens = parse_formula(match.group(1))
        mul = match.group(2)
        if not mul:
            mul = 1
        multiply(parens, int(mul))
        formDict.update(parens)
    # REST
    for match in re.finditer(r"(\(?)([A-Z][a-z]?)(\d*)(\)?)", formula):
        left, elem, mul, right = match.groups()
        if left or right:
            continue
        if not mul:
            mul = 1
        if elem in formDict:
            formDict[elem] += int(mul)
        else:
            formDict[elem] = int(mul)

    return formDict


def formula_to_string(formDict):
    s = ""
    for key, value in formDict.items():
        if value == 1:
            s += key
        elif value > 1:
            s += f"{key}{value}"
    return s


def substract_formulae(minuend, substrahend):
    result = copy.deepcopy(minuend)  # we make a deepcopy to not alter the minuend
    for key, value in substrahend.items():
        result[key] -= value
    return result


def string_formula_substraction(minuend, substrahend):
    return formula_to_string(
        substract_formulae(parse_formula(minuend), parse_formula(substrahend))
    )


"""A dictionary to look up protecting group properties"""
pg_dict = {
    "boc": ("C5H8O2", 100.0524),
    "cbz": ("C8H6O2", 134.0368),
    "tbu": ("C4H8", 56.0626),
    "tms": ("C3H8Si", 72.0395),
}
