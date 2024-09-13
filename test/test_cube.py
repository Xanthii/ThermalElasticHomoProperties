import os
import sys
sys.path.append(r'D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\src')


from lattice_homo.unit_cell import BCCUnitCell, CubeUnitCell
from lattice_homo.problem import ElasticThermalStrategy

geo_params1={'dimX': 1.0}
material_properties = {
    "name": "Material-1",
    "type": "ISOTROPIC",
    "E": 1.0,
    "nu": 0.3,
    "CTE": 5.0e-6,
    "k": 100
}


mesh_params = {
    "deviationFactor":  0.1,
    "minSizeFactor" : 0.1,
    "meshSens" : 1e-7,
    "meshSize": 0.2
}


model_name='Cube_comparision'
cube = CubeUnitCell(model_name, geo_params1, mesh_params, material_properties)
hp = ElasticThermalStrategy()
hp.solve(cube)


# hp.post_process(bcc, display_in_viewport=True)
