# the method of modeling and meshing OCT Unit Cell with 
# 1/8 reduced-scaled model was proved to be correct
# include homogenization analyst process.
from abaqus import *
from abaqusConstants import *
import part
import sketch
import math
import numpy 
import mesh
# Geometry params global
dimX, dimY, dimZ = 1.0, 1.0, 1.0 

inner_radius_xy = 0.085
inner_radius_yz = 0.085
inner_radius_zx = 0.085

outer_radius_xy = 0.085
outer_radius_yz = 0.085
outer_radius_zx = 0.085


# half of the length of edges
xDimension = dimX * 0.5
yDimension = dimY * 0.5
zDimension = dimZ * 0.5
modelName = 'OCT_Cell'
Mdb()
mdb.Model(name=modelName)
if(modelName!='Model-1'):
    del mdb.models['Model-1']

cellVolume = dimX*dimY*dimZ
model = mdb.models[modelName]
# create sketches
skCylinderInner_xy = model.ConstrainedSketch(
    name='CylinderInner_xy',
    sheetSize=inner_radius_xy*2)
skCylinderInner_yz = model.ConstrainedSketch(
    name='CylinderInner_yz',
    sheetSize=inner_radius_yz*2)
skCylinderInner_zx = model.ConstrainedSketch(
    name='CylinderInner_zx',
    sheetSize=inner_radius_zx*2)

skCylinderExter_xy = model.ConstrainedSketch(
    name='CylinderExter_xy',
    sheetSize=outer_radius_xy*2)
skCylinderExter_yz = model.ConstrainedSketch(
    name='CylinderExter_yz',
    sheetSize=outer_radius_yz*2)
skCylinderExter_zx = model.ConstrainedSketch(
    name='CylinderExter_zx',
    sheetSize=outer_radius_zx*2)

# the Rectangular Prism cavity
skCuboid_1 = model.ConstrainedSketch(
    name='Cuboid_1',
    sheetSize=max(xDimension, yDimension))
# the Rectangular Prism cutter
skCuboid_2 = model.ConstrainedSketch(
    name='Cuboid_2',
    sheetSize=5.0*max(xDimension, yDimension))

# create inner cylinder parallel to the xy plane.
skCylinderInner_xy.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(inner_radius_xy, 0.0))
height = (xDimension**2 + yDimension**2)**0.5
partCylinderInner_xy = model.Part(
    name='CylinderInner_xy',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderInner_xy.BaseSolidExtrude(sketch=skCylinderInner_xy, depth=height)
# create inner cylinder parallel to the yz plane.
skCylinderInner_yz.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(inner_radius_yz, 0.0))
height = (zDimension**2 + yDimension**2)**0.5
partCylinderInner_yz = model.Part(
    name='CylinderInner_yz',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderInner_yz.BaseSolidExtrude(sketch=skCylinderInner_yz, depth=height)
# create inner cylinder parallel to the zx plane.
skCylinderInner_zx.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(inner_radius_zx, 0.0))
height = (xDimension**2 + zDimension**2)**0.5
partCylinderInner_zx = model.Part(
    name='CylinderInner_zx',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderInner_zx.BaseSolidExtrude(sketch=skCylinderInner_zx, depth=height)


# create exter cylinder parallel to the xy plane.
skCylinderExter_xy.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_xy, 0.0))
height = (xDimension**2 + yDimension**2)**0.5
partCylinderExter_xy = model.Part(
    name='CylinderExter_xy',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderExter_xy.BaseSolidExtrude(sketch=skCylinderExter_xy, depth=height)
# create exter cylinder parallel to the yz plane.
skCylinderExter_yz.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_yz, 0.0))
height = (xDimension**2 + yDimension**2)**0.5
partCylinderExter_yz = model.Part(
    name='CylinderExter_yz',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderExter_yz.BaseSolidExtrude(sketch=skCylinderExter_yz, depth=height)
# create exter cylinder parallel to the zx plane.
skCylinderExter_zx.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(outer_radius_zx, 0.0))
height = (xDimension**2 + yDimension**2)**0.5
partCylinderExter_zx = model.Part(
    name='CylinderExter_zx',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderExter_zx.BaseSolidExtrude(sketch=skCylinderExter_zx, depth=height)

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

# exter struts
instCylinderExter_xy = assembly.Instance(name='CylinderExter_xy',
    part=partCylinderExter_xy, dependent=OFF)
target_vector = numpy.array([xDimension, yDimension, 0])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderExter_xy',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderExter_xy',),
                vector=(0.0, 0.0, zDimension))

instCylinderExter_yz = assembly.Instance(name='CylinderExter_yz',
    part=partCylinderExter_yz, dependent=OFF)
target_vector = numpy.array([0, yDimension, zDimension])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderExter_yz',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderExter_yz',),
                vector=(xDimension, 0.0, 0.0))

instCylinderExter_zx = assembly.Instance(name='CylinderExter_zx',
    part=partCylinderExter_zx, dependent=OFF)
target_vector = numpy.array([xDimension, 0, zDimension])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderExter_zx',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderExter_zx',),
                vector=(0.0, yDimension, 0.0))


# inner struts
instCylinderInner_xy = assembly.Instance(name='CylinderInner_xy',
    part=partCylinderInner_xy, dependent=OFF)
target_vector = numpy.array([-xDimension, yDimension, 0])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderInner_xy',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderInner_xy',),
                vector=(xDimension, 0.0, 0.0))

instCylinderInner_yz = assembly.Instance(name='CylinderInner_yz',
    part=partCylinderInner_yz, dependent=OFF)
target_vector = numpy.array([0, -yDimension, zDimension])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                (numpy.linalg.norm(initial_vector) 
                * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderInner_yz',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderInner_yz',),
                vector=(0.0, yDimension, 0.0))

instCylinderInner_zx = assembly.Instance(name='CylinderInner_zx',
    part=partCylinderInner_zx, dependent=OFF)
target_vector = numpy.array([xDimension, 0, -zDimension])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderInner_zx',), 
                axisPoint=(0.0, 0.0, 0.0), 
                axisDirection=(rotation_axis[0], rotation_axis[1], rotation_axis[2]), 
                angle=rotation_angle)
assembly.translate(instanceList=('CylinderInner_zx',),
                vector=(0.0, 0.0, zDimension))
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
inList.append(instCylinderInner_xy)
inList.append(instCylinderInner_yz)
inList.append(instCylinderInner_zx)
inList.append(instCylinderExter_xy)
inList.append(instCylinderExter_yz)
inList.append(instCylinderExter_zx)

partCylinders = assembly.PartFromBooleanMerge(
    name='Cylinders',
    instances=inList,
    keepIntersections=False,
    domain=GEOMETRY)

instCylinders = assembly.Instance(name='Cylinders',
    part=partCylinders, dependent=OFF)

del inList
del assembly.instances['CylinderExter_xy']
del assembly.instances['CylinderExter_yz']
del assembly.instances['CylinderExter_zx']
del assembly.instances['CylinderInner_xy']
del assembly.instances['CylinderInner_yz']
del assembly.instances['CylinderInner_zx']

del model.parts['CylinderExter_xy']
del model.parts['CylinderExter_yz']
del model.parts['CylinderExter_zx']
del model.parts['CylinderInner_xy']
del model.parts['CylinderInner_yz']
del model.parts['CylinderInner_zx']

# Cut the cylinders
# create cutter
partCutter = assembly.PartFromBooleanCut(
    name='Cutter',
    instanceToBeCut=instCuboid_2,
    cuttingInstances=(instCuboid_1,), originalInstances=DELETE)
instCutter = assembly.Instance(name='Cutter',
    part=partCutter, dependent=OFF)
# cut
partOCT_Cell_1_8 = assembly.PartFromBooleanCut(
    name='OCT_Cell_1_8',
    instanceToBeCut=instCylinders,
    cuttingInstances=(instCutter,), originalInstances=DELETE)

# instOCT_Cell_1_8 = assembly.Instance(name='OCT_Cell_1_8',
#     part=partOCT_Cell_1_8, dependent=ON)

del model.parts['Cuboid_1']
del model.parts['Cuboid_2']
del model.parts['Cutter']
del model.parts['Cylinders']

# Meshing
deviationFactor = 0.1
minSizeFactor = 0.1
meshSize = xDimension / 12

partOCT_Cell_1_8.seedPart(size=meshSize, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor) 
partOCT_Cell_1_8.setMeshControls(elemShape=TET, regions=
        partOCT_Cell_1_8.cells, technique=FREE)

elemType1 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
partOCT_Cell_1_8.setElementType(regions=(partOCT_Cell_1_8.cells,), elemTypes=(elemType1,))
partOCT_Cell_1_8.generateMesh()


partMeshOCT_1_8= partOCT_Cell_1_8.PartFromMesh(name='MeshPartOCT_1_8', copySets=True)
# Mirror operation in the z direction
partMeshOCT_1_4 = model.PartFromMeshMirror(name='MeshPartOCT_1_4',part=partMeshOCT_1_8,
point1=(xDimension, 0, 0),
point2=(xDimension, 0, 1))
# Mirror operation in the y direction
partMeshOCT_1_2 = model.PartFromMeshMirror(name='MeshPartOCT_1_2',part=partMeshOCT_1_4,
point1=(xDimension, 0, 0),
point2=(xDimension, -1, 0))
# Mirror operation in the x direction
partMeshOCT = model.PartFromMeshMirror(name='MeshPartOCT',part=partMeshOCT_1_2,
point1=(0, yDimension, 0),
point2=(-1, yDimension, 0))

OCTInst = assembly.Instance(name='OCT_Cell', part=partMeshOCT, dependent=ON)

del model.parts['OCT_Cell_1_8']
del model.parts['MeshPartOCT_1_8']
del model.parts['MeshPartOCT_1_4']
del model.parts['MeshPartOCT_1_2']

instanceName = 'OCT_Cell'