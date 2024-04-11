import warnings
from typing import List, Union, Tuple, Sequence

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdChemReactions import (
    ChemicalReaction,
    ReactionFromSmarts,
    ReactionToSmiles,
)

from src.util.rdkit_util import (
    map_reactions,
)


class SFReactionGenerator:
    # reactants to product A
    _rxn_abt = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:1](=[O:2])-[N:3]-[C:4]-[C:5]-[c:6]1:[n:7]:[c:8]2:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:2:[s:14]:1"
    _rxn_th = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[c:7]1:[n:9]:[n:10]:[c:6](-[C:5]-[C:4]-[N:3]-[C:1]=[O:2]):[s:8]:1"

    # reactants to product B
    _rxn_abt_B = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:1](=[O:2])-[N:3]-[C:4]-[C:5]-[C:6](-C(=O)-[OH1])-1-[NH1:7]-[c:8]2:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:2-[S:14]-1"
    _rxn_th_B = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[C:7]-1=[N:9]-[N:10]-[C:6](-C(=O)-[OH1])(-[C:5]-[C:4]-[N:3]-[C:1]=[O:2])-[S:8]-1"

    # reactants to product C
    _rxn_abt_C = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:1]=1-[N:3]-[C:4]-[C:5]-[C:6](-C(=O)-[O-])-2-[N+:7]=1-[c:8]:3:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:3-[S:14]-2"
    _rxn_th_C = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[C:7]-1=[N:9]-[N+:10]=2-[C:6](-C(=O)-[O-])(-[C:5]-[C:4]-[N:3]-[C:1]=2)-[S:8]-1"

    # reactants to product D
    _rxn_abt_D = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[c:1]:1:[n:7]:[c:8]:2:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:2:[s:14]:1"
    _rxn_th_D = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[c:1]:1:[n:10]:[n:9]:[c:7]:[s:8]:1"

    # reactants to product E
    _rxn_abt_E = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[NH2:7]-[c:8]:1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[S:14]-[S:14]-[c:13]:1:[c:12]:[c:11]:[c:10]:[c:9]:[c:8]:1-[NH2:7]"
    _rxn_th_E = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[c:7]:1:[n:9]:[n:9]:[c:7]:[s:8]:1"

    # reactants to product F
    _rxn_abt_F = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:1](=[O:2])-[N:3]-[C:4]-[C:5]-[C:6](=O)-C(=O)-[OH1]"
    _rxn_th_F = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[OH1]-C(=O)-[C:6](=O)-[C:5]-[C:4]-[N:3]-[C:1]=[O:2]"

    # reactants to product G
    _rxn_abt_G = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:1](=[O:2])-[N:3]-[C:4]-[C:5]-[C:6](=O)-[OH1]"
    _rxn_th_G = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5]-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[OH1]-[C:6](=O)-[C:5]-[C:4]-[N:3]-[C:1]=[O:2]"

    # reactants to product H
    _rxn_abt_H = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5](-[H])-[C:4]-[NH1:3]-O-1.[NH2:7]-[c:8]1:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:1-[SH1:14]>>[C:4]=[C:5]-[c:6]:1:[n:7]:[c:8]2:[c:9]:[c:10]:[c:11]:[c:12]:[c:13]:2:[s:14]:1"
    _rxn_th_H = "[$(B(-F)(-F)-F)]-[$(C-[#6])X3:1]=[O:2].O=C1-O-[$(C2CCCCC2)]-O-[C:6]-1-1-[C:5](-[H])-[C:4]-[NH1:3]-O-1.[C:7](=[S:8])-[NH1:9]-[NH2:10]>>[s:8]:1:[c:7]:[n:9]:[n:10]:[c:6]:1-[C:5]=[C:4]"

    # product A to reactants
    _backwards_rxn_abt = "[$(C-[#6]):1](=O)-[NR0]-[C:2]-[C:3]-c1nc2[c:4][c:5][c:6][c:7]c2s1>>F-[B-](-F)(-F)-[C:1]=O.O=C1-O-C2(-C-C-C-C-C-2)-O-C-1-1-[C:3]-[C:2]-N-O1.N-c1:[c:4]:[c:5]:[c:6]:[c:7]:c:1-S"
    _backwards_rxn_th = "[c:4]1nnc(-[C:3]-[C:2]-[NR0]-[$(C-[#6]):1]=O)s1>>F-[B-](-F)(-F)-[C:1]=O.O=C1-O-C2(-C-C-C-C-C-2)-O-C-1-1-[C:3]-[C:2]-N-O1.[C:4](=S)-N-N"

    def __init__(self):
        # Note: the reaction Sanitization warnings are due to the fact that the reaction templates contain dummy atoms
        # that are not mapped to the products (b/c they are not contained in the products ffs)
        # I don't get why RDKit does this as it is an intended and necessary use-case, but as usual their documentation
        # does not help with finding the reason.
        # Sample warnings:
        # Could not find RLabel mapping for atom: 0 in template: 0
        # Mismatched potential rlabels: 1 unmapped reactant dummy atom rlabels,
        # 0 unmappped (sic!) product dummy atom rlabels

        self.backwards_reactions = {
            "abt": ReactionFromSmarts(self._backwards_rxn_abt),
            "th": ReactionFromSmarts(self._backwards_rxn_th),
        }

        self.forward_reactions = {
            "abt": {
                "A": ReactionFromSmarts(self._rxn_abt),
                "B": ReactionFromSmarts(self._rxn_abt_B),
                "C": ReactionFromSmarts(self._rxn_abt_C),
                "D": ReactionFromSmarts(self._rxn_abt_D),
                "E": ReactionFromSmarts(self._rxn_abt_E),
                "F": ReactionFromSmarts(self._rxn_abt_F),
                "G": ReactionFromSmarts(self._rxn_abt_G),
                "H": ReactionFromSmarts(self._rxn_abt_H),
            },
            "th": {
                "A": ReactionFromSmarts(self._rxn_th),
                "B": ReactionFromSmarts(self._rxn_th_B),
                "C": ReactionFromSmarts(self._rxn_th_C),
                "D": ReactionFromSmarts(self._rxn_th_D),
                "E": ReactionFromSmarts(self._rxn_th_E),
                "F": ReactionFromSmarts(self._rxn_th_F),
                "G": ReactionFromSmarts(self._rxn_th_G),
                "H": ReactionFromSmarts(self._rxn_th_H),
            },
        }

        # n.b. we do not use SanitizeRxn here because for the reaction to form product H, we need to specify the loss
        # of a hydrogen atom from the monomer. SanitizeRxn removes the explicit hydrogen atom (WHY GREG, WHY??).
        for rxn in self.backwards_reactions.values():
            rxn.Initialize()
            # SanitizeRxn(rxn)
        for rxn in self.forward_reactions["abt"].values():
            rxn.Initialize()
            # SanitizeRxn(rxn)
        for rxn in self.forward_reactions["th"].values():
            rxn.Initialize()
            # SanitizeRxn(rxn)

    def _try_reaction(
        self, product_type: str, product_mol: Chem.Mol
    ) -> Union[Tuple[List[List[Chem.Mol]], str], Tuple[None, None]]:
        reactants = self.backwards_reactions[product_type].RunReactants((product_mol,))
        if len(reactants) > 0:
            # sanitize before returning
            [[Chem.SanitizeMol(m) for m in pair] for pair in reactants]
            return reactants, product_type
        return None, None

    def generate_reactants(self, product: Union[str, Chem.Mol]) -> List[Chem.Mol]:
        """
        Takes an SF reaction product (product A) and generates the reactants leading to this product.
        Note that the monomer will have an unspecified chiral center at the spiro carbon atom of the keto-acid
        protecting group. This chiral center is not transferred to the product.

        Args:
            product (str or Chem.Mol): Product of a Synthetic Fermentation reaction.

        Returns:
            list: A tuple of 3 reactants (as Mol objects) in the order initiator, monomer, terminator.

        Raises:
            ValueError: If the given product is not a valid Synthetic Fermentation reaction product.
            RuntimeError: If more than one possible reactions is found. This is a sanity check and should never happen.
                We include it because SMARTS pattern matching can sometimes give unforeseen results due to symmetry.
        """
        if isinstance(product, str):
            product_mol = Chem.MolFromSmiles(product)
        else:
            product_mol = Chem.Mol(product)

        # we try applying both reaction templates. We stop as soon as a template works
        # abt
        reactants, product_type = self._try_reaction("abt", product_mol)
        if not reactants:
            # th
            reactants, product_type = self._try_reaction("th", product_mol)
        if not reactants:
            raise ValueError(
                "The given product is not a valid Synthetic Fermentation reaction product."
            )

        # sanity check: there should never be more than one possible reaction
        if len(reactants) > 1:
            raise RuntimeError("More than one possible reaction found.")

        return reactants[0]

    def generate_product(
        self, reactants: List[Chem.Mol], product_type: str = "A"
    ) -> Chem.Mol:
        """
        Generates the product of an SF reaction given the reactants.

        Args:
            reactants (list): Reactants of the Synthetic Fermentation reaction.
                Expects an initiator, monomer, and terminator in this order.
            product_type (str): Type of the product to generate. Can be any from "A" to "H". Defaults to "A".

        Returns:
            Chem.Mol: Product of the SF reaction.
        """

        if product_type == "H":
            monomer = rdmolops.AddHs(reactants[1])
        else:
            monomer = reactants[1]
        initiator = reactants[0]
        terminator = reactants[2]
        try:
            prods = self.forward_reactions["abt"][product_type].RunReactants(
                [initiator, monomer, terminator]
            )
        except KeyError:
            raise ValueError(
                f"Invalid product type '{product_type}'. Must be between 'A' – 'H'."
            )
        if len(prods) == 0:
            try:
                prods = self.forward_reactions["th"][product_type].RunReactants(
                    [initiator, monomer, terminator]
                )
            except KeyError:
                raise ValueError(
                    f"Invalid product type '{product_type}'. Must be between 'A' – 'H'."
                )

        if len(prods) == 0:
            if product_type == "H":
                warnings.warn(
                    "No product H found. This can be expected behavior if the input monomer is a beta-2 spiro monomer"
                )
                return None
            else:
                raise RuntimeError("No product found.")

        elif len(prods) > 1:
            if (
                product_type == "H" and len(prods) == 2
            ):  # H may give two (identical) products if the monomer has two Hs the beta-2 carbon
                pass
            else:
                raise RuntimeError("More than one product found.")

        prod = prods[0][0]
        if product_type == "H":
            prod = rdmolops.RemoveHs(prod)

        Chem.SanitizeMol(prod)
        return prod

    def generate_reaction(
        self,
        reactants: Sequence[Chem.Mol],
    ) -> Union[ChemicalReaction, Tuple[ChemicalReaction, Chem.Mol]]:
        """
        Generates an atom-mapped, unbalanced reaction from the reactants and product type.
        Note that the monomer reactant will have an unspecified chiral center at the spiro carbon atom of the keto-acid.

        Args:
            reactants (Sequence): Reactants of the Synthetic Fermentation reaction.
                Expects an initiator, monomer, and terminator in this order.

        Returns:
            ChemicalReaction: A ChemicalReaction object representing the SF reaction.

        Raises:
            ValueError: If not exactly one reaction is found for the given reactants.
        """

        try:
            reaction = map_reactions(
                self.forward_reactions["abt"]["A"], (reactants,), error_level="error"
            )[
                0
            ]  # we expect exactly one reaction (or an exception)
        except RuntimeError as e1:
            # if this didn't work, try TH reaction template
            try:
                reaction = map_reactions(
                    self.forward_reactions["th"]["A"], (reactants,), error_level="error"
                )[
                    0
                ]  # we expect exactly one reaction (or an exception)
            except RuntimeError as e2:
                raise RuntimeError(
                    f"Encountered Errors on ABT and TH reaction templates while generating reaction for building blocks: "
                    f"'{Chem.MolToSmiles(reactants[0])}', "
                    f"'{Chem.MolToSmiles(reactants[1])}', "
                    f"and '{Chem.MolToSmiles(reactants[2])}'.\n"
                    f"Original error messages: {e1}\n{e2}"
                )

        return reaction

    def get_reaction_smiles(
        self,
        product: Union[str, Chem.Mol],
    ) -> str:
        """
        Generates the atom-mapped SF reactionSMILES leading to the given product.
        Note that the monomer reactant will have an unspecified chiral center at the spiro carbon atom of the keto-acid.
        This means that the reactionSMILES generated here may not be identical to the reactionSMILES in our
        "curated data", which was generated by a different method.
        Differences may also occur in the exact numbers used for atom-mapping, but not in the mapping itself.
        For ML, this is inconsequential, as the CGR generally invariant to the mapping numbers and, using our "custom"
        featurizer also invariant to chirality.

        Args:
            product (str or Chem.Mol): Product (type A) of an SF reaction.

        Returns:
            str: reactionSMILES leading to the given product.
        """
        reactants = self.generate_reactants(product)

        try:
            reaction = self.generate_reaction(reactants)

        except (RuntimeError, ValueError) as e:
            raise RuntimeError(
                f"Encountered Error while generating reaction for product '{product}'.\n"
                f"Original error message: {e}"
            )

        return ReactionToSmiles(reaction)
