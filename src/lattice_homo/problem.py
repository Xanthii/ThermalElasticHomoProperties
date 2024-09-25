from abc import ABC, abstractmethod
from lattice_homo.unit_cell import UnitCell
from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.abaqus_factory import AbaqusUtilityFactory
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
    def solve(cls, model_name: str, unit_cell: UnitCell):
        pass

    @classmethod
    @abstractmethod
    def _process_load_vec(cls, load_vec, vol_unitcell):
        pass

    @classmethod
    def reset(cls):
        '''
        Reset the following class variables.
        '''
        cls._nodes_matched_cache = None
        cls._is_other_analysis_done = False
        AbaqusUtilities.reset_abaqus_environment()


    @classmethod
    def _common_solve_operations(cls, unit_cell:UnitCell, load_vec, suffix, display_in_viewport=True):
        """Handle common solving operations."""
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
        load_vec = cls._process_load_vec(load_vec, vol_unitcell)
        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)

        # the first time to analyze this model.
        # the following operations are performed only once.
        # Default operation
        if not cls._is_other_analysis_done:
            abaqus_utilities.create_model(model_name, geo_params_scaled, mesh_params_scaled)
            abaqus_utilities.assign_material(model_name, material_properties, scale)
            abaqus_utilities.reference_points(model_name)
            cls._nodes_matched_cache = abaqus_utilities.setup_nodes(model_name, mesh_sens, dims)
        # the other analysis has been done.
        # Previous settings need to be cleared.
        else: 
            abaqus_utilities.clear_model_analysis_setup(model_name)
        
        nodes_matched = cls._nodes_matched_cache

        # dynamically 
        getattr(abaqus_utilities, f"set_meshType_{suffix}")(model_name, eleType)
        getattr(abaqus_utilities, f"set_PBC_{suffix}")(model_name, nodes_matched, dims)
        getattr(abaqus_utilities, f"create_steps_{suffix}")(model_name)
        getattr(abaqus_utilities, f"apply_loads_{suffix}")(model_name, load_vec)
        result_odb = getattr(abaqus_utilities, f"submit_job_{suffix}")(model_name, load_vec, display_in_viewport=display_in_viewport)
        effective_props = getattr(abaqus_utilities, f"calc_effective_props_{suffix}")(vol_unitcell, load_vec, result_odb, scale)
        getattr(abaqus_utilities, f"save_effective_props_{suffix}")(model_name, effective_props)
        getattr(abaqus_utilities, f"create_report_json_{suffix}")(model_name, unit_cell, effective_props, scale)

        cls._is_other_analysis_done = True


class ThermalElasticStrategy(ProblemTypeStrategy):
    @classmethod
    def _process_load_vec(cls, load_vec, vol_unitcell):
        # specially handle load_vec
        load_vec = [i * vol_unitcell for i in load_vec]
        load_vec[-1] = load_vec[-1] / vol_unitcell
        return load_vec


    @classmethod
    def solve(cls, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
        suffix = 'te'
        cls._common_solve_operations(unit_cell, load_vec, suffix, display_in_viewport)
       

    @classmethod
    def post_process(cls, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
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
    def _process_load_vec(cls, load_vec, vol_unitcell):
        # specially handle load_vec
        return [i * vol_unitcell for i in load_vec]

    @classmethod
    def solve(cls, unit_cell: UnitCell, load_vec=[1.0]*3, display_in_viewport=True):
        suffix = 'hc'
        cls._common_solve_operations(unit_cell, load_vec, suffix, display_in_viewport)
        
    
    @classmethod
    def post_process(cls, unit_cell: UnitCell, load_vec=[1.0]*3, display_in_viewport=True):
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

        abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)
        
        
        result_odb = abaqus_utilities.get_odb_hc(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props_hc(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props_hc(model_name, effective_props)
        abaqus_utilities.create_report_json_hc(model_name, unit_cell, effective_props,  scale)


