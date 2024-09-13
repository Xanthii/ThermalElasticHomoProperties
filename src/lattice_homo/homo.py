from abaqus import *
from abc import ABC, abstractmethod
from abaqus_utilities import AbaqusUtilities_BCC

import numpy as np
class HomogenizationStrategy(ABC):
    def __init__(self, unit_cells):
        """
        初始化均匀化问题

        :param unit_cells: 一个包含多个单胞实例的列表
        """
        self.unit_cells = unit_cells


    def _compute_scale(self, min_size: float):
        geomScale = np.log10(min_size)
        scale = 10.0 ** int(1 - geomScale)
        return scale

        
    @abstractmethod
    def compute_effective_properties(self):
        pass


class ElasticThermalHomogenizationStrategy(HomogenizationStrategy):
    def __init__(self, unit_cells):
        """
        初始化均匀化问题
        """
        super().__init__(unit_cells)
        self.model_name = self.unit_cells.model_name
        

        # compute geometry scale
        min_size = min(self.unit_cells.get_geo_params().values())
        self.scale = self._compute_scale(min_size)

        self.geo_params_scaled = {k: v * self.scale for k, v in self.unit_cells.get_geo_params().items()}

        # the meshing params will be packed in a dict or dataclass

        self.mesh_params_scaled = {
            k: (v * self.scale if k in ["meshSens", "meshSize"] else v)
            for k, v in self.unit_cells.get_mesh_params().items()
        }
        
     

        self.fea_params ={}
  




    def step(self):
        # create unitcell geometry model and generate meshing
        AbaqusUtilities_BCC.create_BCC_unitcells(self.model_name, self.geo_params_scaled, self.mesh_params_scaled)
        AbaqusUtilities_BCC.material_definition(self.model_name, self.unit_cells.get_material_properties(), self.scale)
        x_dim, y_dim, z_dim = self.geo_params_scaled['dimX'], self.geo_params_scaled['dimY'], self.geo_params_scaled['dimZ']
        AbaqusUtilities_BCC.mesh_handling(self.model_name, self.mesh_params_scaled['meshSens'], x_dim, y_dim, z_dim )
        AbaqusUtilities_BCC.reference_points(self.model_name)
        AbaqusUtilities_BCC.set_PBC(model_name, nodes_matched, x_dim, y_dim, z_dim)

        # AbaqusUtilities_BCC.notes_matching(self.model_name)
        # AbaqusUtilities_BCC.reference_points(self.model_name)
        return nodes_matched
        


    def compute_effective_properties(self):
            pass




