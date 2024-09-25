from lattice_homo.oct_abaqus_utilities import OCTAbaqusUtilities
from lattice_homo.bcc_abaqus_utilities import BCCAbaqusUtilities
from lattice_homo.cuboid_abaqus_utilities import CuboidAbaqusUtilities
from lattice_homo.iso_abaqus_utilities import ISOAbaqusUtilities
from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.unit_cell import UnitCell, BCCUnitCell, CuboidUnitCell, OCTUnitCell, ISOUnitCell

class AbaqusUtilityFactory:
    @staticmethod
    def get_utilities(unit_cell: UnitCell):
        if isinstance(unit_cell, BCCUnitCell):
            return BCCAbaqusUtilities()
        if isinstance(unit_cell, CuboidUnitCell):
            return CuboidAbaqusUtilities()
        if isinstance(unit_cell, ISOUnitCell):
            return ISOAbaqusUtilities()
        elif isinstance(unit_cell, OCTUnitCell):
            return OCTAbaqusUtilities()
        else:
            raise AbaqusUtilities()

