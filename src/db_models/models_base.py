"""
Base models to interact with the database.

To extend functionality, subclass these models in models_custom.py
"""

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import declarative_base, relationship
from rdkit import Chem

Base = declarative_base()


class BuildingBlockBase(Base):
    """
    A building block has few necessary _input_ properties:
        - name (unique)
        - SMILES
        - an nickname for the context of the HTLibrary (unique, can be user-defined or auto-generated)
        - a reactant class (can be "generic")
    A few things will always be generated:
        - id (unique, generated in DB)
        - molblock representation of the building block
        - bottles: Reference to bottles containing this

    """

    __tablename__ = "buildingblock"

    id = Column(Integer, primary_key=True)
    name = Column(String)
    smiles = Column(String)
    nickname = Column(String)
    reactant_class = Column(String)
    molblock = Column(String)
    customization = Column(String)

    bottles = relationship("BottleBase", back_populates="buildingblock")

    __mapper_args__ = {
        "polymorphic_on": customization,
        "polymorphic_identity": "buildingblockbase",
    }

    def __init__(
        self,
        name,
        smiles,
        nickname="auto",
        reactant_class="generic",
        id_int=None,
    ):
        self.name = name
        self.smiles = smiles
        self.nickname = nickname
        self.reactant_class = reactant_class
        self.id = id_int
        self.mol = self.__get_mol()

    def __repr__(self):
        return f"BuildingBlock(id={self.id}, name={self.name}, smiles={self.smiles})"

    def __get_mol(self):
        """Return a molblock from SMILES"""
        mol = Chem.MolFromSmiles(self.smiles)
        Chem.SanitizeMol(mol)
        return Chem.MolToMolBlock(mol)


class BottleBase(Base):
    """
    Base class for chemical bottles.
    A bottle needs at least:
    - id (unique primary key)
    - buildingblock_id (ForeignKey referencing BuildingBlock.id)
    - barcode (nullable)
    -
    """

    __tablename__ = "bottle"

    id = Column(Integer, primary_key=True)
    buildingblock_id = Column(Integer, ForeignKey("buildingblock.id"))
    barcode = Column(Integer)
    customization = Column(String)

    buildingblock = relationship("BuildingBlockBase", back_populates="bottles")
    # TODO are these relationships a problem? Do I need to overwrite them in child classes?

    __mapper_args__ = {
        "polymorphic_on": customization,
        "polymorphic_identity": "bottlebase",
    }

    def __repr__(self):
        return f"Bottle(id={self.id}, barcode={self.barcode}, buildingblock={self.buildingblock})"
