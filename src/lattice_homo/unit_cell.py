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
        - "inner_radius_diag" : float
        - "outer_radius_x" : float
        - "outer_radius_y" : float
        - "outer_radius_z" : float
    '''
    def __post_init__(self):
        super().__post_init__()
        required_keys = ["dimX", "inner_radius_diag", "outer_radius_x"]
        for key in required_keys:
            if key not in self.geo_params:
                raise ValueError(f"Missing required geo_param: {key}")

        if 'dimY' not in self.geo_params:
            self.geo_params['dimY'] = self.geo_params['dimX']
        if 'dimZ' not in self.geo_params:
            self.geo_params['dimZ'] = self.geo_params['dimX']

        if 'outer_radius_y' not in self.geo_params:
            self.geo_params['outer_radius_y'] = self.geo_params['outer_radius_x']
        if 'outer_radius_z' not in self.geo_params:
            self.geo_params['outer_radius_z'] = self.geo_params['outer_radius_x']

        if 'name' not in self.material_properties:
            self.material_properties['name'] = 'Material-1'

        if "num_node_shortest_side" in self.mesh_params.keys():
            num_node_shortest_side = self.mesh_params['num_node_shortest_side'] 
            self.mesh_params['meshSize'] = min(self.geo_params['dimX'],self.geo_params['dimY'],self.geo_params['dimZ']) \
                                            / num_node_shortest_side
        self._apply_geo_scaling()

@dataclass
class ISOUnitCell(BCCUnitCell):
    '''
    A subclass of UnitCell representing a ISO unit cell.
    This class implements specific geometric and mesh parameter handling
    tailored to ISO unit cell structures.

    geo_params : dict 
        A dictionary of geometric parameters with the following structure:
        - "dimX" : float
        - "dimY" : float
        - "dimZ" : float
        - "inner_radius_diag" : float
        - "outer_radius_x" : float
        - "outer_radius_y" : float
        - "outer_radius_z" : float
        - "inner_radius_x" : float
        - "inner_radius_y" : float
        - "inner_radius_z" : float
    '''
    def __post_init__(self):        
        required_keys = ["inner_radius_x"]
        for key in required_keys:
            if key not in self.geo_params:
                raise ValueError(f"Missing required geo_param: {key}")
        if 'inner_radius_y' not in self.geo_params:
            self.geo_params['inner_radius_y'] = self.geo_params['inner_radius_x']
        if 'inner_radius_z' not in self.geo_params:
            self.geo_params['inner_radius_z'] = self.geo_params['inner_radius_x']
        super().__post_init__()
        
@dataclass       
class OCTUnitCell(UnitCell):
    '''
    A subclass of UnitCell representing a Octet truss unit cell.
    This class implements specific geometric and mesh parameter handling
    tailored to OCT unit cell structures.

    geo_params : dict 
        A dictionary of geometric parameters with the following structure:
        - "dimX" : float
        - "dimY" : float
        - "dimZ" : float
        - "outer_radius_xy" : float
        - "outer_radius_yz" : float
        - "outer_radius_zx" : float
        - "inner_radius_xy" : float
        - "inner_radius_yz" : float
        - "inner_radius_zx" : float
        
        "outer_radius_xy" means the outer strut parallel to the xy plane.
        "inner_radius_xy" means the inner strut parallel to the xy plane.


    '''
    def __post_init__(self):        
        required_keys = ["dimX", "inner_radius_x"]
        for key in required_keys:
            if key not in self.geo_params:
                raise ValueError(f"Missing required geo_param: {key}")

        if 'dimY' not in self.geo_params:
            self.geo_params['dimY'] = self.geo_params['dimX']
        if 'dimZ' not in self.geo_params:
            self.geo_params['dimZ'] = self.geo_params['dimX']

        if 'inner_radius_y' not in self.geo_params:
            self.geo_params['inner_radius_y'] = self.geo_params['inner_radius_x']
        if 'inner_radius_z' not in self.geo_params:
            self.geo_params['inner_radius_z'] = self.geo_params['inner_radius_x']
        super().__post_init__()

@dataclass
class CuboidUnitCell(UnitCell):
    # need to be extended to outer_radius variables.
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