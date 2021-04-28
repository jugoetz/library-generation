import copy
import re


def parse_formula(formula: str) -> dict:  # Formula Parsing by Aditya Matam
    def multiply(formula: dict, mul: int) -> None:
        for key in formula: formula[key] *= mul

    formDict = {}
    # PARENS
    for match in re.finditer(r"\((.*?)\)(\d*)", formula):
        parens = parse_formula(match.group(1))
        mul = match.group(2)
        if not mul: mul = 1
        multiply(parens, int(mul))
        formDict.update(parens)
    # REST
    for match in re.finditer(r"(\(?)([A-Z][a-z]?)(\d*)(\)?)", formula):
        left, elem, mul, right = match.groups()
        if left or right: continue
        if not mul: mul = 1
        if elem in formDict:
            formDict[elem] += int(mul)
        else:
            formDict[elem] = int(mul)

    return formDict


def formula_to_string(formDict):
    s = ''
    for key, value in formDict.items():
        if value == 1:
            s += key
        elif value > 1:
            s += f'{key}{value}'
    return s


def substract_formulae(minuend, substrahend):
    result = copy.deepcopy(minuend)  # we make a deepcopy to not alter the minuend
    for key, value in substrahend.items():
        result[key] -= value
    return result


def string_formula_substraction(minuend, substrahend):
    return formula_to_string(substract_formulae(parse_formula(minuend), parse_formula(substrahend)))
