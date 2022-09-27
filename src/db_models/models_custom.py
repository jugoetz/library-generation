"""
Custom models to interact with the database.
"""

from sqlalchemy import Column, String, Integer, Date, ForeignKey, Float
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from src.db_models.models_base import Base, BottleBase, BuildingBlockBase


class MyBuildingBlock(BuildingBlockBase):
    """
    Subclass of BuildingBlock with extended attributes.
    Additional attributes given by user input:
        - comment
        - category
    Additional attributes generated from mol:
        - weigh-in
        - protecting groups TODO
        - lcms masses TODO
        - lcms formulae TODO

    Also, inputs can be subject to filtering rules. TODO
    """
    __mapper_args__ = {
        'polymorphic_identity': 'mybuildingblock'
    }

    category = Column(String)
    comment = Column(String)
    weigh_in = Column(String)

    def __init__(self, name, smiles, category, nickname='auto', reactant_class='generic', id_int=None, comment=None,
                 ):
        super().__init__(name, smiles, nickname, reactant_class, id_int
                         )
        self.category = category
        self.comment = comment
        self.weigh_in = self.__get_weigh_in()

    def __repr__(self):
        return f'BuildingBlock(id={self.id}, name={self.name}, smiles={self.smiles}, category={self.category}, weigh-in={self.weigh_in})'

    def __get_weigh_in(self):
        """Calculate weigh-in for preparing 0.1 mL of a 0.05 M solution from MolWt"""
        weight = MolWt(Chem.MolFromMolBlock(self.mol))
        return round(weight * 1e-4 * 0.05 * 1000, 2)


class MyBottle(BottleBase):
    """
    Custom class derived off the base class for chemical bottles
    Additional info:
    - mass / mass unit
    """
    __mapper_args__ = {
        'polymorphic_identity': 'mybottle'
    }

    mass = Column(Float)
    mass_unit = Column(String)

    def __repr__(self):
        return f'Bottle(id={self.id}, barcode={self.barcode}, buildingblock={self.buildingblock})'
