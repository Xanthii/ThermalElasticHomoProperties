from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Any
import numpy as np


@dataclass
class UnitCell:
    '''
    A base class representing the basic properties and methods of a unit cell
    in multi-scale structural design. The class defines parameters for geometry, 
    mesh configuration, material properties, and a scaling mechanism.

    Parameters
    ----------
        model_name : str
            The name of the cell
        geo_params : dict 
            A dictionary of geometric parameters according to cell type
        mesh_params : dict 
            A dictionary of mesh parameters with the following structure:
            - "deviationFactor" : float
            - "minSizeFactor" : float
            - "meshSens" : float
            - "meshSize" : float
            - "eleType" : str
        material_properties : dict 
            A dictionary of material properties with the following structure:
            - "name" : str
            - "type" : str
            - "E" : float
            - "nu" : float
            - "CTE" : float
            - "k" : float
    '''
 
    model_name: str
    geo_params: Dict[str, float]
    mesh_params: Dict[str, Any]
    material_properties: Dict[str, Any]

    def __post_init__(self):
        # Merge provided material_properties with defaults
        default_material_properties = {
            "name": "Material-1",
            "type": "ISOTROPIC",
            "E": 1.0,
            "nu": 0.3,
            "CTE": 0.0,
            "k": 0.0
        }
        self.material_properties = {**default_material_properties, **self.material_properties}


    def get_geo_params(self) -> Dict[str, Dict[str, Any]]:
        return self.geo_params

    def get_material_properties(self) -> Dict[str, Dict[str, Any]]:
        return self.material_properties
    
    def get_mesh_params(self) -> Dict[str, Dict[str, Any]]:
        return self.mesh_params

    def get_scale(self) -> float:
        return self.scale

    def _apply_geo_scaling(self):
        min_size = min(self.geo_params.values())
        self.scale = UnitCell._compute_scale(min_size)
        # self.geo_scaled = {k: v * self.scale for k, v in self.geo_params.items()}

    @staticmethod
    def _compute_scale(min_size: float):
        geomScale = np.log10(min_size)
        scale = 10.0 ** int(1 - geomScale)
        return scale
    
    
@dataclass
class BCCUnitCell(UnitCell):
    '''
    A subclass of UnitCell representing a Body-Centered Cubic (BCC) unit cell.
    This class implements specific geometric and mesh parameter handling
    tailored to BCC unit cell structures.

    geo_params : dict 
        A dictionary of geometric parameters with the following structure:
        - "dimX" : float
        - "dimY" : float
        - "dimZ" : float
        - "radiusInter" : float
        - "radiusExter_x" : float
        - "radiusExter_y" : float
        - "radiusExter_z" : float
    '''
    def __post_init__(self):
        super().__post_init__()

        required_keys = ["dimX", "radiusInter", "radiusExter_x"]
        for key in required_keys:
            if key not in self.geo_params:
                raise ValueError(f"Missing required geo_param: {key}")

        if 'dimY' not in self.geo_params:
            self.geo_params['dimY'] = self.geo_params['dimX']
        if 'dimZ' not in self.geo_params:
            self.geo_params['dimZ'] = self.geo_params['dimX']

        if 'radiusExter_y' not in self.geo_params:
            self.geo_params['radiusExter_y'] = self.geo_params['radiusExter_x']
        if 'radiusExter_z' not in self.geo_params:
            self.geo_params['radiusExter_z'] = self.geo_params['radiusExter_x']

        if 'name' not in self.material_properties:
            self.material_properties['name'] = 'Material-1'

        if "num_node_shortest_side" in self.mesh_params.keys():
            num_node_shortest_side = self.mesh_params['num_node_shortest_side'] 
            self.mesh_params['meshSize'] = min(self.geo_params['dimX'],self.geo_params['dimY'],self.geo_params['dimZ']) \
                                            / num_node_shortest_side
        self._apply_geo_scaling()


class FCCUnitCell(UnitCell):
    pass
        
class OCTUnitCell(UnitCell):
    pass


@dataclass
class CuboidUnitCell(UnitCell):
    '''
    A subclass of UnitCell representing a Cuboid unit cell.
    This class implements specific geometric and mesh parameter handling
    tailored to Cuboid unit cell structures.

    geo_params : dict 
        A dictionary of geometric parameters with the following structure:
        - "dimX" : float
        - "dimY" : float
        - "dimZ" : float
    '''
    def __post_init__(self):
        required_keys = ["dimX"]
        for key in required_keys:
            if key not in self.geo_params:
                raise ValueError(f"Missing required geo_param: {key}")
        if 'dimY' not in self.geo_params:
            self.geo_params['dimY'] = self.geo_params['dimX']
        if 'dimZ' not in self.geo_params:
            self.geo_params['dimZ'] = self.geo_params['dimX']

        if 'name' not in self.material_properties:
            self.material_properties['name'] = 'Material-1'

        if "num_node_shortest_side" in self.mesh_params.keys():
            num_node_shortest_side = self.mesh_params['num_node_shortest_side'] 
            self.mesh_params['meshSize'] = min(self.geo_params['dimX'],self.geo_params['dimY'],self.geo_params['dimZ']) \
                                            / num_node_shortest_side

        self._apply_geo_scaling()