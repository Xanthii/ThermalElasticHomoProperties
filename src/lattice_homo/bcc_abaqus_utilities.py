from abaqus import *
from abaqusConstants import *
import mesh
import numpy as np
from typing import Dict, Any

from lattice_homo.abaqus_utilities import AbaqusUtilities
from lattice_homo.abaqus_utilities import  ThermalElasticAbaqusUtilities
from lattice_homo.abaqus_utilities import  HeatConductionAbaqusUtilities

class BCCAbaqusUtilities(ThermalElasticAbaqusUtilities, HeatConductionAbaqusUtilities):
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
        inner_radius_diag = geo_params_scaled['inner_radius_diag']
        outer_radius_x, outer_radius_y, outer_radius_z= (geo_params_scaled[k] for k in ['outer_radius_x', 'outer_radius_y', 'outer_radius_z'])
        eleType = AbaqusUtilities._get_eleType(mesh_params_scaled['eleType'])
        Mdb()
        model = mdb.Model(name=model_name)
        if(model_name != 'Model-1'):
            del mdb.models['Model-1']

        # create sketches
        skCylinderInner_diag = model.ConstrainedSketch(
            name='CylinderInner_diag',
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

        # the Rectangular Prism cavity
        skCuboid_1 = model.ConstrainedSketch(
            name='Cuboid_1',
            sheetSize=max(xDimension, yDimension))

        # the Rectangular Prism cutter
        skCuboid_2 = model.ConstrainedSketch(
            name='Cuboid_2',
            sheetSize=5.0*max(xDimension, yDimension))

        # create inter cylinder
        skCylinderInner_diag.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(inner_radius_diag, 0.0))
        height = (xDimension**2 + yDimension**2 + zDimension**2)**0.5
        partCylinderInner_diag = model.Part(
            name='CylinderInner_diag',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCylinderInner_diag.BaseSolidExtrude(sketch=skCylinderInner_diag, depth=height)

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

        # create Rectangular Prism cavity
        skCuboid_1.rectangle(
            point1=(0.0, 0.0),
            point2=(xDimension, yDimension))
        partCuboid_1 = model.Part(
            name='Cuboid_1',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCuboid_1.BaseSolidExtrude(sketch=skCuboid_1, depth=zDimension)

        # create Rectangular Prism cutter
        skCuboid_2.rectangle(
            point1=(0.0, 0.0),
            point2=(5*xDimension, 5*yDimension))
        partCuboid_2 = model.Part(
            name='Cuboid_2',
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        partCuboid_2.BaseSolidExtrude(sketch=skCuboid_2, depth=5*zDimension)

        # Instantiate the parts
        # Then rotate and translate
        assembly = model.rootAssembly
        # inter struts
        instCylinderInner_diag = assembly.Instance(name='CylinderInner_diag',
            part=partCylinderInner_diag, dependent=OFF)
        # calculate the rotation angle and axis

        target_vector = np.array([xDimension, yDimension, zDimension])
        initial_vector = np.array([0, 0, zDimension]) 
        rotation_axis = np.cross(initial_vector, target_vector)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        rotation_angle = degrees(acos(np.dot(initial_vector, target_vector) / 
                            (np.linalg.norm(initial_vector) 
                            * np.linalg.norm(target_vector))))
        assembly.rotate(instanceList=('CylinderInner_diag',), 
                        axisPoint=(0.0, 0.0, 0.0), 
                        axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                        angle=rotation_angle)

        # exter struts
        instCylinderExter_x = assembly.Instance(name='CylinderExter_x',
            part=partCylinderExter_x, dependent=OFF)
        assembly.rotate(instanceList=('CylinderExter_x',), 
                        axisPoint=(0.0, 0.0, 0.0), 
                        axisDirection=(0.0, 1.0, 0.0), 
                        angle=90.0) 
        assembly.translate(instanceList=('CylinderExter_x',),
                        vector=(0.0, yDimension, zDimension))

        instCylinderExter_y = assembly.Instance(name='CylinderExter_y',
            part=partCylinderExter_y, dependent=OFF)
        assembly.rotate(instanceList=('CylinderExter_y',), 
                        axisPoint=(0.0, 0.0, 0.0), 
                        axisDirection=(1.0, 0.0, 0.0), 
                        angle=-90.0)
        assembly.translate(instanceList=('CylinderExter_y',),
                        vector=(xDimension, 0.0, zDimension))

        instCylinderExter_z = assembly.Instance(name='CylinderExter_z',
            part=partCylinderExter_z, dependent=OFF)
        assembly.translate(instanceList=('CylinderExter_z',),
                        vector=(xDimension, yDimension, 0.0))

        # cavity
        instCuboid_1 = assembly.Instance(name='Cuboid_1',
            part=partCuboid_1, dependent=OFF)
        # cutter box
        instCuboid_2 = assembly.Instance(name='Cuboid_2',
            part=partCuboid_2, dependent=OFF)
        assembly.translate(instanceList=('Cuboid_2',),
                        vector=(-2*xDimension, -2*yDimension, -2*zDimension))

        # Merge the four cylinder components.
        inList=[]
        inList.append(instCylinderInner_diag)
        inList.append(instCylinderExter_x)
        inList.append(instCylinderExter_y)
        inList.append(instCylinderExter_z)
        partCylinders = assembly.PartFromBooleanMerge(
            name='Cylinders',
            instances=inList,
            keepIntersections=False,
            domain=GEOMETRY)

        instCylinders = assembly.Instance(name='Cylinders',
            part=partCylinders, dependent=OFF)

        del inList
        del assembly.instances['CylinderInner_diag']
        del assembly.instances['CylinderExter_x']
        del assembly.instances['CylinderExter_y']
        del assembly.instances['CylinderExter_z']

        del model.parts['CylinderInner_diag']
        del model.parts['CylinderExter_x']
        del model.parts['CylinderExter_y']
        del model.parts['CylinderExter_z']

        # Cut the cylinders
        # create cutter
        partCutter = assembly.PartFromBooleanCut(
            name='Cutter',
            instanceToBeCut=instCuboid_2,
            cuttingInstances=(instCuboid_1,), originalInstances=DELETE)
        instCutter = assembly.Instance(name='Cutter',
            part=partCutter, dependent=OFF)
        # cut
        partBCC_Cell_1_8 = assembly.PartFromBooleanCut(
            name='BCC_Cell_1_8',
            instanceToBeCut=instCylinders,
            cuttingInstances=(instCutter,), originalInstances=DELETE)

        # instBCC_Cell_1_8 = assembly.Instance(name='BCC_Cell_1_8',
        #     part=partBCC_Cell_1_8, dependent=ON)

        del model.parts['Cuboid_1']
        del model.parts['Cuboid_2']
        del model.parts['Cutter']
        del model.parts['Cylinders']


        # ______________________________ #
        #      meshing 
        # ______________________________ #
        deviationFactor = mesh_params_scaled['deviationFactor']
        minSizeFactor = mesh_params_scaled['minSizeFactor']
        meshSize = mesh_params_scaled['meshSize']
      
        partBCC_Cell_1_8.seedPart(size=meshSize, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor) 
        partBCC_Cell_1_8.setMeshControls(elemShape=TET, regions=
                partBCC_Cell_1_8.cells, technique=FREE)

        elemType1 = mesh.ElemType(elemCode=eleType, elemLibrary=STANDARD)
        partBCC_Cell_1_8.setElementType(regions=(partBCC_Cell_1_8.cells,), elemTypes=(elemType1,))
        partBCC_Cell_1_8.generateMesh()


        partMeshBCC_1_8= partBCC_Cell_1_8.PartFromMesh(name='MeshPartBCC_1_8', copySets=True)

        partMeshBCC_1_4 = model.PartFromMeshMirror(name='MeshPartBCC_1_4',part=partMeshBCC_1_8,
        point1=(xDimension-outer_radius_z,yDimension-outer_radius_z,0),
        point2=(xDimension-outer_radius_z,yDimension-outer_radius_z,1))

        partMeshBCC_1_2 = model.PartFromMeshMirror(name='MeshPartBCC_1_2',part=partMeshBCC_1_4,
        point1=(xDimension-outer_radius_y,0,zDimension-outer_radius_y),
        point2=(xDimension-outer_radius_y,-1,zDimension-outer_radius_y))


        part_name = AbaqusUtilities._get_part_name(model_name)
        partMeshBCC = model.PartFromMeshMirror(name=part_name,part=partMeshBCC_1_2,
        point1=(0,yDimension-outer_radius_x,zDimension-outer_radius_x),
        point2=(-1,yDimension-outer_radius_x,zDimension-outer_radius_x))

        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        BCCInst = assembly.Instance(name=instance_name, part=partMeshBCC, dependent=ON)

        del model.parts['BCC_Cell_1_8']
        del model.parts['MeshPartBCC_1_8']
        del model.parts['MeshPartBCC_1_4']
        del model.parts['MeshPartBCC_1_2']
