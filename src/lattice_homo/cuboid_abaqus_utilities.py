from abaqus import *
from abaqusConstants import *
import mesh
import numpy as np
from typing import Dict, Any

from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.abaqus_utilities import  ThermalElasticAbaqusUtilities
from lattice_homo.abaqus_utilities import  HeatConductionAbaqusUtilities

class CuboidAbaqusUtilities(ThermalElasticAbaqusUtilities, HeatConductionAbaqusUtilities):
    @staticmethod
    def create_model(model_name: str, geo_params_scaled: Dict[str, float], mesh_params_scaled: Dict[str, float]):
        """
        Abaqus-python script to create a meshed instance named 'BCC_Cell'.
        Assume that the origin is at the center of the cell.

        geo_params_scaled: Contains the same keys as in unitcell, but the values() were scaled.
        Assigning material section is done here.
        """
        # half of the original length
        xDimension, yDimension, zDimension = (geo_params_scaled[k] * 0.5 for k in ['dimX', 'dimY', 'dimZ'])

        Mdb()
        model = mdb.Model(name=model_name)
        if(model_name != 'Model-1'):
            del mdb.models['Model-1']

        # create sketches
        skCube_1_8 = model.ConstrainedSketch(
            name='Cube_1_8',
            sheetSize=max(xDimension, yDimension))

        skCube_1_8.rectangle(
            point1=(0.0, 0.0),
            point2=(xDimension, yDimension))
        partCube_1_8 = model.Part(
            name='Cube_Cell_1_8',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCube_1_8.BaseSolidExtrude(sketch=skCube_1_8, depth=zDimension)
        assembly = model.rootAssembly
      # ______________________________ #
        deviationFactor = mesh_params_scaled['deviationFactor']
        minSizeFactor = mesh_params_scaled['minSizeFactor']
        meshSize = mesh_params_scaled['meshSize']
      
        partCube_1_8.seedPart(size=meshSize, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor) 
        partCube_1_8.setMeshControls(elemShape=TET, regions=
                partCube_1_8.cells, technique=FREE)

        elemType1 = mesh.ElemType(elemCode=eleType, elemLibrary=STANDARD)
        partCube_1_8.setElementType(regions=(partCube_1_8.cells,), elemTypes=(elemType1,))
        partCube_1_8.generateMesh()

        partMeshCube_1_8= partCube_1_8.PartFromMesh(name='MeshPartCube_1_8', copySets=True)

        partMeshCube_1_4 = model.PartFromMeshMirror(name='MeshPartCube_1_4',part=partMeshCube_1_8,
        point1=(xDimension * 0.5, yDimension * 0.5, 0),
        point2=(xDimension * 0.5, yDimension * 0.5,-1))

        partMeshCube_1_2 = model.PartFromMeshMirror(name='MeshPartCube_1_2',part=partMeshCube_1_4,
        point1=(xDimension * 0.5, 0, zDimension * 0.5),
        point2=(xDimension * 0.5,-1, zDimension * 0.5))


        part_name = AbaqusUtilities._get_part_name(model_name)
        partMeshCube = model.PartFromMeshMirror(name=part_name,part=partMeshCube_1_2,
        point1=(0, yDimension * 0.5, zDimension * 0.5),
        point2=(-1, yDimension * 0.5, zDimension * 0.5))

        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        CubeInst = assembly.Instance(name=instance_name, part=partMeshCube, dependent=ON)

        del model.parts['Cube_Cell_1_8']
        del model.parts['MeshPartCube_1_8']
        del model.parts['MeshPartCube_1_4']
        del model.parts['MeshPartCube_1_2']
