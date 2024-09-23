from abc import ABC, abstractmethod
from lattice_homo.unit_cell import UnitCell
from lattice_homo.abaqus_utilities import AbaqusUtilityFactory
import math
class ProblemTypeStrategy(ABC):
    '''
    provide the unified interface.
    '''

    # Class variables, storing intermediate calculation results
    _nodes_matched_cache = None
    _is_other_analysis_done = False

    @classmethod
    @abstractmethod
    def solve(self, model_name: str, unit_cell: UnitCell):
        pass

class ElasticThermalStrategy(ProblemTypeStrategy):
    @classmethod
    def solve(cls, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
        # Abaqus script dealing with elastic thermodynamics
        model_name = unit_cell.model_name
        scale = unit_cell.get_scale()
        geo_params_scaled = {k: v * scale for k, v in unit_cell.get_geo_params().items()}
        
        mesh_params_scaled = {
            k: (v * scale if k in ["meshSens", "meshSize"] else v)
            for k, v in unit_cell.get_mesh_params().items()
            }
        material_properties =  unit_cell.get_material_properties()
        mesh_sens = mesh_params_scaled['meshSens']
        
        dims = [geo_params_scaled[k] for k in ['dimX', 'dimY', 'dimZ']]
        vol_unitcell = math.prod(dims)
        load_vec = [i * vol_unitcell for i in load_vec]
        load_vec[-1] = load_vec[-1] / vol_unitcell

        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)
        
        if not cls._is_other_analysis_done:
            abaqus_utilities.create_model(model_name, geo_params_scaled, mesh_params_scaled)
            abaqus_utilities.assign_material(model_name, material_properties, scale)
            abaqus_utilities.reference_points(model_name)
            cls._nodes_matched_cache = abaqus_utilities.setup_nodes(model_name, mesh_sens, dims)

        nodes_matched = cls._nodes_matched_cache
        abaqus_utilities.set_PBC_et(model_name, nodes_matched, dims)
        abaqus_utilities.create_steps_et(model_name)
        abaqus_utilities.apply_loads_et(model_name, load_vec)
        result_odb = abaqus_utilities.submit_job_et(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props_et(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props_et(model_name, effective_props)
        abaqus_utilities.create_report_json_et(model_name, unit_cell, effective_props,  scale)



    def post_process(self, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
        # Abaqus script dealing with elastic thermodynamics
        model_name = unit_cell.model_name
        scale = unit_cell.get_scale()
        geo_params_scaled = {k: v * scale for k, v in unit_cell.get_geo_params().items()}
        
        mesh_params_scaled = {
            k: (v * scale if k in ["meshSens", "meshSize"] else v)
            for k, v in unit_cell.get_mesh_params().items()
            }
  
        dims = [geo_params_scaled[k] for k in ['dimX', 'dimY', 'dimZ']]
        vol_unitcell = math.prod(dims)
        load_vec = [i * vol_unitcell for i in load_vec]
        load_vec[-1] = load_vec[-1] / vol_unitcell

        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)
        
        
        result_odb = abaqus_utilities.get_odb(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props(model_name, effective_props)
        abaqus_utilities.create_report_json(model_name, unit_cell, effective_props, scale)




class HeatConductionStrategy(ProblemTypeStrategy):
    @classmethod
    def solve(cls, unit_cell: UnitCell, load_vec=[1.0]*3, display_in_viewport=True):
        model_name = unit_cell.model_name
        scale = unit_cell.get_scale()
        geo_params_scaled = {k: v * scale for k, v in unit_cell.get_geo_params().items()}
        
        mesh_params_scaled = {
            k: (v * scale if k in ["meshSens", "meshSize"] else v)
            for k, v in unit_cell.get_mesh_params().items()
            }
        material_properties =  unit_cell.get_material_properties()
        mesh_sens = mesh_params_scaled['meshSens']
        eleType = mesh_params_scaled['eleType']
        dims = [geo_params_scaled[k] for k in ['dimX', 'dimY', 'dimZ']]
        vol_unitcell = math.prod(dims)
        load_vec = [i * vol_unitcell for i in load_vec]

        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)

        
        if not cls._is_other_analysis_done:
            abaqus_utilities.create_model(model_name, geo_params_scaled, mesh_params_scaled)
            abaqus_utilities.assign_material(model_name, material_properties, scale)
            abaqus_utilities.reference_points(model_name)
            cls._nodes_matched_cache = abaqus_utilities.setup_nodes(model_name, mesh_sens, dims)
        print('here')
        abaqus_utilities._set_meshType_hc(model_name, eleType)
        nodes_matched = cls._nodes_matched_cache
        abaqus_utilities.set_PBC_hc(model_name, nodes_matched, dims)
        abaqus_utilities.create_steps_hc(model_name)
        abaqus_utilities.apply_loads_hc(model_name, load_vec)
        result_odb = abaqus_utilities.submit_job_hc(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props_hc(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props_hc(model_name, effective_props)
        abaqus_utilities.create_report_json_hc(model_name, unit_cell, effective_props,  scale)
    

    def post_process(self, unit_cell: UnitCell, load_vec=[1.0]*3, display_in_viewport=True):
        # Abaqus script dealing with elastic thermodynamics
        model_name = unit_cell.model_name
        scale = unit_cell.get_scale()
        geo_params_scaled = {k: v * scale for k, v in unit_cell.get_geo_params().items()}
        
        mesh_params_scaled = {
            k: (v * scale if k in ["meshSens", "meshSize"] else v)
            for k, v in unit_cell.get_mesh_params().items()
            }
        material_properties =  unit_cell.get_material_properties()
        mesh_sens = mesh_params_scaled['meshSens']
        dims = [geo_params_scaled[k] for k in ['dimX', 'dimY', 'dimZ']]
        vol_unitcell = math.prod(dims)
        load_vec = [i * vol_unitcell for i in load_vec]

        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)
        
        
        result_odb = abaqus_utilities.get_odb_hc(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props_hc(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props_hc(model_name, effective_props)
        abaqus_utilities.create_report_json_hc(model_name, unit_cell, effective_props, scale)
