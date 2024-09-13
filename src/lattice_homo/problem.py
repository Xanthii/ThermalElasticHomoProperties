from abc import ABC, abstractmethod
from lattice_homo.unit_cell import UnitCell
from lattice_homo.abaqus_utilities import AbaqusUtilityFactory
import math
class ProblemTypeStrategy(ABC):
    '''
    provide the unified interface.
    '''
    @abstractmethod
    def solve(self, model_name: str, unit_cell: UnitCell):
        pass

class ElasticThermalStrategy(ProblemTypeStrategy):
    def solve(self, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
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
        # include geometry model and mesh model
        
        
        abaqus_utilities.create_model(model_name, geo_params_scaled, mesh_params_scaled)
        abaqus_utilities.assign_material(model_name, material_properties, scale)
        abaqus_utilities.reference_points(model_name)
        abaqus_utilities.setup_mesh_pbc(model_name, mesh_sens, dims)
        abaqus_utilities.create_steps(model_name)
        abaqus_utilities.apply_loads(model_name, load_vec)
        result_odb = abaqus_utilities.submit_job(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props(model_name, effective_props)
        abaqus_utilities.create_report_json(model_name, unit_cell, effective_props,  scale)



    def post_process(self, unit_cell: UnitCell, load_vec=[1.0]*7, display_in_viewport=True):
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
        
        
        result_odb = abaqus_utilities.get_odb(model_name, load_vec, display_in_viewport)
        effective_props = abaqus_utilities.calc_effective_props(vol_unitcell, load_vec, result_odb, scale)
        abaqus_utilities.save_effective_props(model_name, effective_props)
        abaqus_utilities.create_report_json(model_name, unit_cell, effective_props, scale)




# class HeatConductionStrategy(ProblemTypeStrategy):
#     def solve(self, model_name: str, unit_cell: UnitCell):
#         # Abaqus script dealing with heat conduction
#         pass
#         abaqus_utilities = AbaqusUtilityFactory.get_utilities(unit_cell)
#         abaqus_utilities.mesh_handling_hc(model_name, ...)
#         abaqus_utilities.set_PBC_hc(model_name, ...)

    