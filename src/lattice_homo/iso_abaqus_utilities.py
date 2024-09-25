from abaqus import *
from abaqusConstants import *
import numpy as np
from typing import Dict, Any

from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.abaqus_utilities import  ThermalElasticAbaqusUtilities
from lattice_homo.abaqus_utilities import  HeatConductionAbaqusUtilities

class ISOAbaqusUtilities(ThermalElasticAbaqusUtilities, HeatConductionAbaqusUtilities):
    @staticmethod
    def create_model(model_name: str, geo_params_scaled: Dict[str, float], mesh_params_scaled: Dict[str, float]):
        xDimension, yDimension, zDimension = (geo_params_scaled[k] * 0.5 for k in ['dimX', 'dimY', 'dimZ'])
        inner_radius_diag = geo_params_scaled['inner_radius_diag']
        outer_radius_x, outer_radius_y, outer_radius_z = (geo_params_scaled[k] for k in ['outer_radius_x', 'outer_radius_y', 'outer_radius_z'])
        inner_radius_x, inner_radius_y, inner_radius_z = (geo_params_scaled[k] for k in ['inner_radius_x', 'inner_radius_y', 'inner_radius_z'])
        eleType = AbaqusUtilities._get_eleType(mesh_params_scaled['eleType'])
        Mdb()
        model = mdb.Model(name=model_name)
        if(model_name != 'Model-1'):
            del mdb.models['Model-1']

        # create sketches
        skCylinderInter_diag = model.ConstrainedSketch(
            name='CylinderInter',
            sheetSize=inner_radius_diag*2)
        skCylinderExter_x = model.ConstrainedSketch(
            name='CylinderExter_x',
            sheetSize=outer_radius_x*2)
        skCylinderExter_y = model.ConstrainedSketch(
            name='CylinderExter_y',
            sheetSize=outer_radius_y*2)
        skCylinderExter_z = model.ConstrainedSketch(
            name='CylinderExter_z',
            sheetSize=outer_radius_z*2)

        skCylinderInner_x = model.ConstrainedSketch(
            name='CylinderInner_x',
            sheetSize=inner_radius_x*2)
        skCylinderInner_y = model.ConstrainedSketch(
            name='CylinderInner_y',
            sheetSize=inner_radius_y*2)
        skCylinderInner_z = model.ConstrainedSketch(
            name='CylinderInner_z',
            sheetSize=inner_radius_z*2)

        # the Rectangular Prism cavity
        skCuboid_1 = model.ConstrainedSketch(
            name='Cuboid_1',
            sheetSize=max(xDimension, yDimension))

        # the Rectangular Prism cutter
        skCuboid_2 = model.ConstrainedSketch(
            name='Cuboid_2',
            sheetSize=5.0*max(xDimension, yDimension))

        # create inter cylinder in diag-direction
        skCylinderInter.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(inner_radius_diag, 0.0))
        height = (xDimension**2 + yDimension**2 + zDimension**2)**0.5
        partCylinderInter = model.Part(
            name='CylinderInter',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCylinderInter.BaseSolidExtrude(sketch=skCylinderInter, depth=height)

        # create exter cylinder in x-direction
        skCylinderExter_x.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_x, 0.0))
        partCylinderExter_x = model.Part(
            name='CylinderExter_x',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCylinderExter_x.BaseSolidExtrude(sketch=skCylinderExter_x, depth=xDimension)
        # create exter cylinder in y-direction
        skCylinderExter_y.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_y, 0.0))
        partCylinderExter_y = model.Part(
            name='CylinderExter_y',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCylinderExter_y.BaseSolidExtrude(sketch=skCylinderExter_y, depth=yDimension)
        # create exter cylinder in z-direction
        skCylinderExter_z.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_z, 0.0))
        partCylinderExter_z = model.Part(
            name='CylinderExter_z',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCylinderExter_z.BaseSolidExtrude(sketch=skCylinderExter_z, depth=zDimension)


