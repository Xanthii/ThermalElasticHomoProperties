from abaqus import *
from abaqusConstants import *
import mesh
import numpy as np
from typing import Dict, Any

from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.abaqus_utilities import  ThermalElasticAbaqusUtilities
from lattice_homo.abaqus_utilities import  HeatConductionAbaqusUtilities

class OCTAbaqusUtilities(ThermalElasticAbaqusUtilities, HeatConductionAbaqusUtilities):
    pass