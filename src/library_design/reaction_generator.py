from typing import List, Union, Tuple, Sequence

from rdkit import Chem
from rdkit.Chem.rdChemReactions import (
    ChemicalReaction,
    ReactionFromSmarts,
    SanitizeRxn,
    ReactionToSmiles,
)

from src.util.rdkit_util import (
    map_reactions,
)


class SFReactionGenerator:
    # reactants to product A
    _rxn_abt = "[$(B(-F)(-F)-F)]-[C:2](-[#6:1])=[O:3].O=C1-O-[$(C2CCCCC2)]-O-[C:7]-1-1-[C:6]-[C:5]-[N:4]-O-1.[N:8]-[c:9]1:[c:10]:[c:11]:[c:12]:[c:13]:[c:14]:1-[S:15]>>[#6:1]-[C:2](=[O:3])-[N:4]-[C:5]-[C:6]-[c:7]1:[n:8]:[c:9]2:[c:10]:[c:11]:[c:12]:[c:13]:[c:14]:2:[s:15]:1."
    _rxn_th = "[$(B(-F)(-F)-F)]-[C:2](-[#6:1])=[O:3].O=C1-O-[$(C2CCCCC2)]-O-[C:7]-1-1-[C:6]-[C:5]-[N:4]-O-1.[#6:8]-[C:9](=[S:10])-[N:11]-[N:12]>>[#6:8]-[c:9]1:[n:11]:[n:12]:[c:7](-[C:6]-[C:5]-[N:4]-[C:2](-[#6:1])=[O:3]):[s:10]:1."

    # product A to reactants
    _backwards_rxn_abt = "[#6:1]-[#6](=O)-[NR0]-[#6:2]-[#6:3]-c1nc2[c:4][c:5][c:6][c:7]c2s1>>F[B-](F)(F)[#6](-[#6:1])=O.O=[#6]1-[#8]C2([#6]-[#6]-[#6]-[#6]-[#6]2)[#8]C11[#6:3]-[#6:2]-[#7]-[#8]1.[#7]-c1[c:4][c:5][c:6][c:7]c1-[#16]"
    _backwards_rxn_th = "[#6:4]-c1nnc(-[C:3]-[C:2]-[NR0]-C(-[#6:1])=O)s1>>F-[B-](-F)(-F)-C(-[#6:1])=O.O=C1-O-C2(-C-C-C-C-C-2)-O-C-1-1-[C:3]-[C:2]-N-O1.[#6:4]-C(=S)-N-N"

    # reactants to product B

    def __init__(self):
        self.backwards_reactions = {
            "abt": ReactionFromSmarts(self._backwards_rxn_abt),
            "th": ReactionFromSmarts(self._backwards_rxn_th),
        }

        self.forward_reactions = {
            "abt": ReactionFromSmarts(self._rxn_abt),
            "th": ReactionFromSmarts(self._rxn_th),
        }

        for rxn in self.backwards_reactions.values():
            rxn.Initialize()
            SanitizeRxn(rxn)
        for rxn in self.forward_reactions.values():
            rxn.Initialize()
            SanitizeRxn(rxn)

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

    def generate_product(self, reactants: Sequence[Chem.Mol]) -> Chem.Mol:
        """
        Generates the product of an SF reaction given the reactants.

        Args:
            reactants (Sequence): Reactants of the Synthetic Fermentation reaction.
                Expects an initiator, monomer, and terminator in this order.

        Returns:
            Chem.Mol: Product of the SF reaction.
        """

        prods = self.forward_reactions["abt"].RunReactants(reactants)
        if len(prods) == 0:
            prods = self.forward_reactions["th"].RunReactants(reactants)

        if len(prods) == 0:
            raise RuntimeError("No product found.")
        elif len(prods) > 1:
            raise RuntimeError("More than one product found.")

        return prods[0][0]

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
                self.forward_reactions["abt"], (reactants,), error_level="error"
            )[
                0
            ]  # we expect exactly one reaction (or an exception)
        except RuntimeError as e1:
            # if this didn't work, try TH reaction template
            try:
                reaction = map_reactions(
                    self.forward_reactions["th"], (reactants,), error_level="error"
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
