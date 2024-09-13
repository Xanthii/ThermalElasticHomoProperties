from unit_cell import BCCUnitCell
from problem import ElasticThermalStrategy

geo_params1={'dimX': 1.0, 
            'radiusInter': 0.05279541,
            'radiusExter_x': 0.07647705
        }
geo_params2={'dimX': 1.0, 
            'radiusInter':  0.14361572265625,
            'radiusExter_x': 0.18463134765625
        }
material_properties = {
    "name": "Material-1",
    "type": "ISOTROPIC",
    "E": 1.0,
    "nu": 0.3,
    "CTE": 0.3
}


mesh_params = {
    "deviationFactor":  0.1,
    "minSizeFactor" : 0.1,
    "meshSens" : 1e-7,
    "meshSize": 0.05
}


model_name='BCC_vertified_2'
bcc = BCCUnitCell(model_name, geo_params1, mesh_params, material_properties)
hp = ElasticThermalStrategy()
# hp.solve(bcc)


hp.post_process(bcc, display_in_viewport=True)
