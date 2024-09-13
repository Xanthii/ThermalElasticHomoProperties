# the method of modeling and meshing BCC Unit Cell with 
# 1/8 reduced-scaled model was proved to be correct
from abaqus import *
from abaqusConstants import *
import part
import sketch
import math
import numpy 
import mesh
# Geometry params global
dimX, dimY, dimZ = 1.0, 1.0, 1.0 
radiusInter = 0.05279541
radiusExter_x = 0.07647705
radiusExter_y = 0.07647705
radiusExter_z = 0.07647705


# half of the length of edges
xDimension = dimX * 0.5
yDimension = dimY * 0.5
zDimension = dimZ * 0.5
modelName = 'BCC_Cell'
Mdb()
mdb.Model(name=modelName)
if(modelName!='Model-1'):
    del mdb.models['Model-1']

cellVolume = dimX*dimY*dimZ
model = mdb.models[modelName]
# create sketches
skCylinderInter = model.ConstrainedSketch(
    name='CylinderInter',
    sheetSize=radiusInter*2)
skCylinderExter_x = model.ConstrainedSketch(
    name='CylinderExter_x',
    sheetSize=radiusExter_x*2)
skCylinderExter_y = model.ConstrainedSketch(
    name='CylinderExter_y',
    sheetSize=radiusExter_y*2)
skCylinderExter_z = model.ConstrainedSketch(
    name='CylinderExter_z',
    sheetSize=radiusExter_z*2)

# the Rectangular Prism cavity
skCuboid_1 = model.ConstrainedSketch(
    name='Cuboid_1',
    sheetSize=max(xDimension, yDimension))

# the Rectangular Prism cutter
skCuboid_2 = model.ConstrainedSketch(
    name='Cuboid_2',
    sheetSize=5.0*max(xDimension, yDimension))

# create inter cylinder
skCylinderInter.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radiusInter, 0.0))
height = (xDimension**2 + yDimension**2 + zDimension**2)**0.5
partCylinderInter = model.Part(
    name='CylinderInter',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderInter.BaseSolidExtrude(sketch=skCylinderInter, depth=height)

# create exter cylinder in x-direction
skCylinderExter_x.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radiusExter_x, 0.0))
partCylinderExter_x = model.Part(
    name='CylinderExter_x',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderExter_x.BaseSolidExtrude(sketch=skCylinderExter_x, depth=xDimension)
# create exter cylinder in y-direction
skCylinderExter_y.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radiusExter_y, 0.0))
partCylinderExter_y = model.Part(
    name='CylinderExter_y',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
partCylinderExter_y.BaseSolidExtrude(sketch=skCylinderExter_y, depth=yDimension)
# create exter cylinder in z-direction
skCylinderExter_z.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radiusExter_z, 0.0))
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
instCylinderInter = assembly.Instance(name='CylinderInter',
    part=partCylinderInter, dependent=OFF)
# calculate the rotation angle and axis

target_vector = numpy.array([xDimension, yDimension, zDimension])
initial_vector = numpy.array([0, 0, zDimension]) 
rotation_axis = numpy.cross(initial_vector, target_vector)
rotation_axis = rotation_axis / numpy.linalg.norm(rotation_axis)
rotation_angle = degrees(acos(numpy.dot(initial_vector, target_vector) / 
                    (numpy.linalg.norm(initial_vector) 
                    * numpy.linalg.norm(target_vector))))
assembly.rotate(instanceList=('CylinderInter',), 
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
inList.append(instCylinderInter)
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
del assembly.instances['CylinderInter']
del assembly.instances['CylinderExter_x']
del assembly.instances['CylinderExter_y']
del assembly.instances['CylinderExter_z']

del model.parts['CylinderInter']
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

# Meshing
deviationFactor = 0.1
minSizeFactor = 0.2
meshSize = xDimension / 6

partBCC_Cell_1_8.seedPart(size=meshSize, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor) 
partBCC_Cell_1_8.setMeshControls(elemShape=TET, regions=
        partBCC_Cell_1_8.cells, technique=FREE)

elemType1 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
partBCC_Cell_1_8.setElementType(regions=(partBCC_Cell_1_8.cells,), elemTypes=(elemType1,))
partBCC_Cell_1_8.generateMesh()


partMeshBCC_1_8= partBCC_Cell_1_8.PartFromMesh(name='MeshPartBCC_1_8', copySets=True)

partMeshBCC_1_4 = model.PartFromMeshMirror(name='MeshPartBCC_1_4',part=partMeshBCC_1_8,
point1=(xDimension-radiusExter_z,yDimension-radiusExter_z,0),
point2=(xDimension-radiusExter_z,yDimension-radiusExter_z,1))

partMeshBCC_1_2 = model.PartFromMeshMirror(name='MeshPartBCC_1_2',part=partMeshBCC_1_4,
point1=(xDimension-radiusExter_y,0,zDimension-radiusExter_y),
point2=(xDimension-radiusExter_y,-1,zDimension-radiusExter_y))

partMeshBCC = model.PartFromMeshMirror(name='MeshPartBCC',part=partMeshBCC_1_2,
point1=(0,yDimension-radiusExter_x,zDimension-radiusExter_x),
point2=(-1,yDimension-radiusExter_x,zDimension-radiusExter_x))

BCCInst = assembly.Instance(name='BCC_Cell', part=partMeshBCC, dependent=ON)

del model.parts['BCC_Cell_1_8']
del model.parts['MeshPartBCC_1_8']
del model.parts['MeshPartBCC_1_4']
del model.parts['MeshPartBCC_1_2']

instanceName = 'BCC_Cell'
# __________________________________________________
# Material
mdb.models[modelName].Material(name="Material-1")

mdb.models[modelName].materials["Material-1"].Elastic(
    type=ISOTROPIC,
    table=((1.0, \
    0.3),))
mdb.models[modelName].materials["Material-1"].Expansion(
                type=ISOTROPIC,
                table=((0.4,),))

mdb.models[modelName].HomogeneousSolidSection(
		name="Material-1",
		material="Material-1",
		thickness=None)	


p = mdb.models[modelName].parts["MeshPartBCC"]
elements = p.elements[:]
region = p.Set(elements=elements, name='ElementSet-1')
p.SectionAssignment(region=region, sectionName="Material-1")



# nodes handling
Nodeset = BCCInst.nodes
meshsens = 1e-6
Max = xDimension
May = yDimension
Maz = zDimension
Mnx = -xDimension
Mny = -yDimension
Mnz = -zDimension
# container 

x=[]
y=[]
z=[]
for i in Nodeset:
    x.append(i.coordinates[0])
    y.append(i.coordinates[1])
    z.append(i.coordinates[2])


vertex1=[]
vertex2=[]
vertex3=[]
vertex4=[]
vertex5=[]
vertex6=[]
vertex7=[]
vertex8=[]

edge_IIIxyz={}
edge_IVxyz={}
edge_IIxyz={}
edge_Ixyz={}
edge_VIIxyz={}
edge_VIIIxyz={}
edge_VIxyz={}
edge_Vxyz={}
edge_XIxyz={}
edge_Xxyz={}
edge_XIIxyz={}
edge_IXxyz={}

faceAxyz={}
faceBxyz={}
faceCxyz={}
faceDxyz={}
faceExyz={}
faceFxyz={}

edge_III=[]
edge_IV=[]
edge_II=[]
edge_I=[]
edge_VII=[]
edge_VI=[]
edge_VIII=[]
edge_V=[]
edge_XI=[]
edge_XII=[]
edge_X=[]
edge_IX=[]

faceA=[]
faceB=[]
faceC=[]
faceD=[]
faceE=[]
faceF=[]

coc1={}
coc2={}
coc3={}
coc4={}
coc5={}
coc6={}
coc7={}
coc8={}
## Identifying boundary nodes ##
for i in Nodeset:
    # interior nodes
    if (Mnx+meshsens) < i.coordinates[0] < (Max-meshsens) and (Mny+meshsens) < i.coordinates[1] < (May-meshsens) and (Mnz+meshsens) < i.coordinates[2] < (Maz-meshsens):
        continue
    # vertex
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
        vertex1.insert(0,i.label)
        coc1[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
        vertex2.insert(0,i.label)
        coc2[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
        vertex3.insert(0,i.label)
        coc3[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
        vertex4.insert(0,i.label)
        coc4[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
        vertex5.insert(0,i.label)
        coc5[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
        vertex6.insert(0,i.label)
        coc6[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
        vertex7.insert(0,i.label)
        coc7[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
        vertex8.insert(0,i.label)
        coc8[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    # edge
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        edge_IIIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        edge_IIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        edge_IVxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        edge_Ixyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
        edge_VIIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
        edge_VIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
        edge_VIIIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
        edge_Vxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        edge_XIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        edge_XIIxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        edge_Xxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        edge_IXxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
    # face
    
    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        faceAxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        faceBxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]  
    if abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        faceCxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
        faceDxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        faceExyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
        faceFxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	

# Matching 
# single-node Set
## Sorting and appending sets ##
# faces
faceA = faceAxyz.keys()
for i in faceA:
    for j in faceBxyz.keys():
        if abs(faceAxyz[i][1] - faceBxyz[j][1])<=meshsens and abs(faceAxyz[i][2] - faceBxyz[j][2])<=meshsens:
            faceB.append(j)

faceC = faceCxyz.keys()
for i in faceC:
    for j in faceDxyz.keys():
        if abs(faceCxyz[i][0] - faceDxyz[j][0]) <=meshsens and abs(faceCxyz[i][2] - faceDxyz[j][2]) <=meshsens:
            faceD.append(j)

faceE = faceExyz.keys()
for i in faceE:
    for j in faceFxyz.keys():
        if abs(faceExyz[i][0] - faceFxyz[j][0])<=meshsens and abs(faceExyz[i][1] - faceFxyz[j][1]) <=meshsens:
            faceF.append(j)


for i,j in zip(faceA,faceB):
    assembly.SetFromNodeLabels(name='Node_FaceA_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_FaceB_%s' % (j), nodeLabels=((instanceName,[j]),))

for i,j in zip(faceC,faceD):
    assembly.SetFromNodeLabels(name='Node_FaceC_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_FaceD_%s' % (j), nodeLabels=((instanceName,[j]),))

for i,j in zip(faceE,faceF):
    assembly.SetFromNodeLabels(name='Node_FaceE_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_FaceF_%s' % (j), nodeLabels=((instanceName,[j]),))

# edges I-IV
edge_I = edge_Ixyz.keys()
for i in edge_I:
    for j in edge_IIxyz.keys():
        if abs(edge_Ixyz[i][2] - edge_IIxyz[j][2])<=meshsens:
            edge_II.append(j)
    for j in edge_IIIxyz.keys():
        if abs(edge_Ixyz[i][2] - edge_IIIxyz[j][2]) <=meshsens:
            edge_III.append(j)
    for j in edge_IVxyz.keys():
        if abs(edge_Ixyz[i][2] - edge_IVxyz[j][2]) <=meshsens:
            edge_IV.append(j)

# edges V-VIII
edge_V = edge_Vxyz.keys()
for i in edge_V:
    for j in edge_VIxyz.keys():
        if abs(edge_Vxyz[i][1] - edge_VIxyz[j][1])<=meshsens:
            edge_VI.append(j)
    for j in edge_VIIxyz.keys():
        if abs(edge_Vxyz[i][1] - edge_VIIxyz[j][1]) <=meshsens:
            edge_VII.append(j)
    for j in edge_VIIIxyz.keys():
        if abs(edge_Vxyz[i][1] - edge_VIIIxyz[j][1]) <=meshsens:
            edge_VIII.append(j)


# edges IX-XII
edge_IX = edge_IXxyz.keys()
for i in edge_IX:
    for j in edge_Xxyz.keys():
        if abs(edge_IXxyz[i][0] - edge_Xxyz[j][0])<=meshsens:
            edge_X.append(j)
    for j in edge_XIxyz.keys():
        if abs(edge_IXxyz[i][0] - edge_XIxyz[j][0]) <=meshsens:
            edge_XI.append(j)
    for j in edge_XIIxyz.keys():
        if abs(edge_IXxyz[i][0] - edge_XIIxyz[j][0]) <=meshsens:
            edge_XII.append(j)

for i,j,k,l in zip(edge_I,edge_II,edge_III,edge_IV):
    assembly.SetFromNodeLabels(name='Node_edge_I_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_edge_II_%s' % (j), nodeLabels=((instanceName,[j]),))
    assembly.SetFromNodeLabels(name='Node_edge_III_%s' % (k), nodeLabels=((instanceName,[k]),))
    assembly.SetFromNodeLabels(name='Node_edge_IV_%s' % (l), nodeLabels=((instanceName,[l]),))

for i,j,k,l in zip(edge_V,edge_VI,edge_VII,edge_VIII):
    assembly.SetFromNodeLabels(name='Node_edge_V_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_edge_VI_%s' % (j), nodeLabels=((instanceName,[j]),))
    assembly.SetFromNodeLabels(name='Node_edge_VII_%s' % (k), nodeLabels=((instanceName,[k]),))
    assembly.SetFromNodeLabels(name='Node_edge_VIII_%s' % (l), nodeLabels=((instanceName,[l]),))

for i,j,k,l in zip(edge_IX,edge_X,edge_XI,edge_XII):
    assembly.SetFromNodeLabels(name='Node_edge_IX_%s' % (i), nodeLabels=((instanceName,[i]),))
    assembly.SetFromNodeLabels(name='Node_edge_X_%s' % (j), nodeLabels=((instanceName,[j]),))
    assembly.SetFromNodeLabels(name='Node_edge_XI_%s' % (k), nodeLabels=((instanceName,[k]),))
    assembly.SetFromNodeLabels(name='Node_edge_XII_%s' % (l), nodeLabels=((instanceName,[l]),))

# vertex
assembly.SetFromNodeLabels(name='Master Node 1', nodeLabels=((instanceName,vertex1),))
assembly.SetFromNodeLabels(name='Master Node 2', nodeLabels=((instanceName,vertex2),))
assembly.SetFromNodeLabels(name='Master Node 3', nodeLabels=((instanceName,vertex3),))
assembly.SetFromNodeLabels(name='Master Node 4', nodeLabels=((instanceName,vertex4),))
assembly.SetFromNodeLabels(name='Master Node 5', nodeLabels=((instanceName,vertex5),))
assembly.SetFromNodeLabels(name='Master Node 6', nodeLabels=((instanceName,vertex6),))
assembly.SetFromNodeLabels(name='Master Node 7', nodeLabels=((instanceName,vertex7),))
assembly.SetFromNodeLabels(name='Master Node 8', nodeLabels=((instanceName,vertex8),))

# from CommonUCFeatures import ReferencePoints
# ReferencePoints(modelName)
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Fy load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Fz load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Shear_yz load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Shear_zx load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Shear_xy load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Tx load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Ty load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))
# Point holding the Tz load.
mdb.models[modelName].rootAssembly.ReferencePoint(
    point=(0.0, 0.0, 0.0))

mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-1', 
    toName='Fx')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-2', 
    toName='Fy')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-3', 
    toName='Fz')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-4', 
    toName='Shear_yz')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-5', 
    toName='Shear_zx')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-6', 
    toName='Shear_xy')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-7', 
    toName='Tx')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-8', 
    toName='Ty')
mdb.models[modelName].rootAssembly.features.changeKey(
    fromName='RP-9', 
    toName='Tz')

keyRefPointList = mdb.models[modelName].rootAssembly.referencePoints.keys()
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[8]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER FX')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[7]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER FY')		
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[6]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER FZ')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[5]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER SHEAR_YZ')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[4]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER SHEAR_ZX')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[3]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='CONSTRAINTS DRIVER SHEAR_XY')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[2]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='Constraints Driver Tx')
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[1]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='Constraints Driver Ty')		
refPoint = (mdb.models[modelName].rootAssembly.referencePoints[keyRefPointList[0]],)
mdb.models[modelName].rootAssembly.Set(
    referencePoints=refPoint,
    name='Constraints Driver Tz')


xDim, yDim, zDim = xDimension, yDimension, zDimension
# PBC
# vertex
# This is the U2-U1=FAB.
mdb.models[modelName].Equation(
    name='U2-U1 Master node on u',
    terms=((1.0,'Master Node 2', 1),
        (-1.0, 'Master Node 1', 1),
        (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U2-U1 Master node on v',
    terms=((1.0,'Master Node 2', 2),
        (-1.0, 'Master Node 1', 2),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U2-U1 Master node on w',
    terms=((1.0,'Master Node 2', 3),
        (-1.0, 'Master Node 1', 3),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))

# This is the U3-U1=FAB+FCD equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U3-U1 Master node on u',
    terms=((1.0,'Master Node 3', 1),
        (-1.0, 'Master Node 1', 1),
        (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U3-U1 Master node on v',
    terms=((1.0,'Master Node 3', 2),
        (-1.0, 'Master Node 1', 2),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
        (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U3-U1 Master node on w',
    terms=((1.0,'Master Node 3', 3),
        (-1.0, 'Master Node 1', 3),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
        (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))	

# This is the U4-U1=FCD equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U4-U1 Master node on u',
    terms=((1.0,'Master Node 4', 1),
        (-1.0, 'Master Node 1', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U4-U1 Master node on v',
    terms=((1.0,'Master Node 4', 2),
        (-1.0, 'Master Node 1', 2),
        (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U4-U1 Master node on w',
    terms=((1.0,'Master Node 4', 3),
        (-1.0, 'Master Node 1', 3),
        (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))

# This is the U5-U1=FEF equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U5-U1 Master node on u',
    terms=((1.0,'Master Node 5', 1),
        (-1.0, 'Master Node 1', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U5-U1 Master node on v',
    terms=((1.0,'Master Node 5', 2),
        (-1.0, 'Master Node 1', 2)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U5-U1 Master node on w',
    terms=((1.0,'Master Node 5', 3),
        (-1.0, 'Master Node 1', 3),
        (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

# This is the U6-U1=FAB+FEF equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U6-U1 Master node on u',
    terms=((1.0,'Master Node 6', 1),
        (-1.0, 'Master Node 1', 1),
        (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U6-U1 Master node on v',
    terms=((1.0,'Master Node 6', 2),
        (-1.0, 'Master Node 1', 2),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U6-U1 Master node on w',
    terms=((1.0,'Master Node 6', 3),
        (-1.0, 'Master Node 1', 3),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
        (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
        
# This is the U7-U1=FAB+FCD+FEF equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U7-U1 Master node on u',
    terms=((1.0,'Master Node 7', 1),
        (-1.0, 'Master Node 1', 1),
        (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U7-U1 Master node on v',
    terms=((1.0,'Master Node 7', 2),
        (-1.0, 'Master Node 1', 2),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
        (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U7-U1 Master node on w',
    terms=((1.0,'Master Node 7', 3),
        (-1.0, 'Master Node 1', 3),
        (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
        (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
        (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))	

# This is the U8-U1=FCD+FEF equation.
# The 'u' equation (Abaqus coordinate 1)
mdb.models[modelName].Equation(
    name='U8-U1 Master node on u',
    terms=((1.0,'Master Node 8', 1),
        (-1.0, 'Master Node 1', 1)))
# The 'v' equation (Abaqus coordinate 2)
mdb.models[modelName].Equation(
    name='U8-U1 Master node on v',
    terms=((1.0,'Master Node 8', 2),
        (-1.0, 'Master Node 1', 2),
        (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
# The 'w' equation (Abaqus coordinate 3)
mdb.models[modelName].Equation(
    name='U8-U1 Master node on w',
    terms=((1.0,'Master Node 8', 3),
        (-1.0, 'Master Node 1', 3),
        (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
        (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

# edges equation
# edge I-IV
for i,j,k,l,count in zip(edge_I,edge_II,edge_III,edge_IV,list(range(1,len(edge_I)+1))):
    nodeSetLabel_I = 'Node_edge_I_%s'%i
    nodeSetLabel_II = 'Node_edge_II_%s'%j
    nodeSetLabel_III = 'Node_edge_III_%s'%k
    nodeSetLabel_IV = 'Node_edge_IV_%s'%l
    # This is the U_II-U_I=FAB equation.
    # The 'u' equation (Abaqus coordinate 1)
    eqEdgeName = 'U_II-U_I edge node '
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_II, 1),
            (-1.0, nodeSetLabel_I, 1),
            (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_II, 2),
            (-1.0, nodeSetLabel_I, 2),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_II, 3),
            (-1.0, nodeSetLabel_I, 3),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))
    # This is the U_III-U_I=FAB+FCD equation.
    eqEdgeName = 'U_III-U_I edge node '
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_III, 1),
            (-1.0, nodeSetLabel_I, 1),
            (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_III, 2),
            (-1.0, nodeSetLabel_I, 2),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
            (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_III, 3),
            (-1.0, nodeSetLabel_I, 3),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
            (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))
    # This is the U_IV-U_I=FCD equation.
    eqEdgeName = 'U_IV-U_I edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_IV, 1),
            (-1.0, nodeSetLabel_I, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_IV, 2),
            (-1.0, nodeSetLabel_I, 2),
            (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_IV, 3),
            (-1.0, nodeSetLabel_I, 3),
            (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))

# edge V-VIII
for i,j,k,l,count in zip(edge_V, edge_VI, edge_VII, edge_VIII, list(range(1,len(edge_V)+1))):
    nodeSetLabel_V = 'Node_edge_V_%s'%i
    nodeSetLabel_VI = 'Node_edge_VI_%s'%j
    nodeSetLabel_VII = 'Node_edge_VII_%s'%k
    nodeSetLabel_VIII = 'Node_edge_VIII_%s'%l
    # This is the U_VI-U_V=FAB equation.
    eqEdgeName = 'U_VI-U_V edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VI, 1),
            (-1.0, nodeSetLabel_V, 1),
            (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VI, 2),
            (-1.0, nodeSetLabel_V, 2),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_VI, 3),
            (-1.0, nodeSetLabel_V, 3),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))
    # This is the U_VII-U_V=FAB+FEF equation.
    eqEdgeName = 'U_VII-U_V edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VII, 1),
            (-1.0, nodeSetLabel_V, 1),
            (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VII, 2),
            (-1.0, nodeSetLabel_V, 2),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_VII, 3),
            (-1.0, nodeSetLabel_V, 3),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
            (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
    # This is the U_VIII-U_V=FEF equation.
    eqEdgeName = 'U_VIII-U_V edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VIII, 1),
            (-1.0, nodeSetLabel_V, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_VIII, 2),
            (-1.0, nodeSetLabel_V, 2)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_VIII, 3),
            (-1.0, nodeSetLabel_V, 3),
            (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

# edge IX-XII
for i,j,k,l,count in zip(edge_IX, edge_X, edge_XI, edge_XII, list(range(1,len(edge_IX)+1))):
    nodeSetLabel_IX = 'Node_edge_IX_%s'%i
    nodeSetLabel_X = 'Node_edge_X_%s'%j
    nodeSetLabel_XI = 'Node_edge_XI_%s'%k
    nodeSetLabel_XII = 'Node_edge_XII_%s'%l
    # This is the U_X-U_IX=FCD equation.
    eqEdgeName = 'U_X-U_IX edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_X, 1),
            (-1.0, nodeSetLabel_IX, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_X, 2),
            (-1.0, nodeSetLabel_IX, 2),
            (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_X, 3),
            (-1.0, nodeSetLabel_IX, 3),
            (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))
    # This is the U_XI-U_IX=FCD+FEF equation.
    eqEdgeName = 'U_XI-U_IX edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_XI, 1),
            (-1.0, nodeSetLabel_IX, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_XI, 2),
            (-1.0, nodeSetLabel_IX, 2),
            (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_XI, 3),
            (-1.0, nodeSetLabel_IX, 3),
            (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
            (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
    # This is the U_XII-U_IX=FEF equation.
    eqEdgeName = 'U_XII-U_IX edge node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqEdgeName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_XII, 1),
            (-1.0, nodeSetLabel_IX, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqEdgeName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_XII, 2),
            (-1.0, nodeSetLabel_IX, 2)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqEdgeName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_XII, 3),
            (-1.0, nodeSetLabel_IX, 3),
            (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

# Equations for the faces.
# This is the U_B-U_A=FAB equation.
for i,j,count in zip(faceA, faceB,list(range(1,len(faceA)+1))):
    nodeSetLabel_A = 'Node_FaceA_%s'%i
    nodeSetLabel_B = 'Node_FaceB_%s'%j
    eqFaceName = 'U_B-U_A face node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqFaceName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_B, 1),
            (-1.0, nodeSetLabel_A, 1),
            (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqFaceName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_B, 2),
            (-1.0, nodeSetLabel_A, 2),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqFaceName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_B, 3),
            (-1.0, nodeSetLabel_A, 3),
            (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))

# This is the U_D-U_C=FCD equation.
for i,j,count in zip(faceC, faceD,list(range(1,len(faceC)+1))):
    nodeSetLabel_C = 'Node_FaceC_%s'%i
    nodeSetLabel_D = 'Node_FaceD_%s'%j
    eqFaceName = 'U_D-U_C face node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqFaceName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_D, 1),
            (-1.0, nodeSetLabel_C, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqFaceName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_D, 2),
            (-1.0, nodeSetLabel_C, 2),
            (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqFaceName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_D, 3),
            (-1.0, nodeSetLabel_C, 3),
            (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))

# This is the U_F-U_E=FEF equation.
for i,j,count in zip(faceE, faceF,list(range(1,len(faceE)+1))):
    nodeSetLabel_E = 'Node_FaceE_%s'%i
    nodeSetLabel_F = 'Node_FaceF_%s'%j
    eqFaceName = 'U_F-U_E face node '
    # The 'u' equation (Abaqus coordinate 1)
    eqName = eqFaceName + str(count) + ' on u'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_F, 1),
            (-1.0, nodeSetLabel_E, 1)))
    # The 'v' equation (Abaqus coordinate 2)
    eqName = eqFaceName + str(count) + ' on v'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0, nodeSetLabel_F, 2),
            (-1.0, nodeSetLabel_E, 2)))
    # The 'w' equation (Abaqus coordinate 3)
    eqName = eqFaceName + str(count) + ' on w'
    mdb.models[modelName].Equation(
        name=eqName,
        terms=((1.0,nodeSetLabel_F, 3),
            (-1.0, nodeSetLabel_E, 3),
            (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

regMasterNode1 = mdb.models[modelName].rootAssembly.sets['Master Node 1']
mdb.models[modelName].DisplacementBC(
    name='Translation stop Master Node 1', 
    createStepName='Initial',
    region=regMasterNode1,
    u1=SET,
    u2=SET,
    u3=SET,
    ur1=UNSET, 
    ur2=UNSET,
    ur3=UNSET,
    amplitude=UNSET,
    distributionType=UNIFORM, 
    localCsys=None)

# print('Define BC Complete')

modelName, previousStep = modelName, 'Initial'
mdb.models[modelName].StaticLinearPerturbationStep(
name='Isothermal linear perturbation step',
description='Computation of effective elastic material properties',
previous=previousStep)

mdb.models[modelName].StaticLinearPerturbationStep(
name='Thermomechanical step',
description='Computation of coefficients of thermal expansion',
previous='Isothermal linear perturbation step')	

# Delete the default history output request.
# History output requests are not compatible with the load cases option.
#mdb.models[modelName].historyOutputRequests['H-Output-1']
mdb.models[modelName].HistoryOutputRequest(name='H-Output-1', 
createStepName='Thermomechanical step', variables=('SJD', 'SJDA', 'SJDT', 
'SJDTA'))

# Stress only field output requests.
mdb.models[modelName].FieldOutputRequest(
name='Standard output request', 
createStepName='Isothermal linear perturbation step',
variables=('S','E','EVOL'))

# Define an an additional field output request for the Fx RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FX']
mdb.models[modelName].FieldOutputRequest(
name='Output Request Fx', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Define an an additional field output request for the Fy RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FY']
mdb.models[modelName].FieldOutputRequest(
name='Ouput Request Fy', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Define an an additional field output request for the Fz RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FZ']
mdb.models[modelName].FieldOutputRequest(
name='Ouput Request Fz', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Define an an additional field output request for the Shear_xy RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_XY']
mdb.models[modelName].FieldOutputRequest(
name='Ouput Request Shear_xy', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Deactivate Shear_xy reporting for the thermomechanical step.
mdb.models[modelName].fieldOutputRequests['Ouput Request Shear_xy'].deactivate('Thermomechanical step')	
# Define an an additional field output request for the Shear_yz RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_YZ']
mdb.models[modelName].FieldOutputRequest(
name='Ouput Request Shear_yz', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Deactivate Shear_yz reporting for the thermomechanical step.
mdb.models[modelName].fieldOutputRequests['Ouput Request Shear_yz'].deactivate('Thermomechanical step')	
# Define an an additional field output request for the Shear_zx RP.
regionDef=mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_ZX']
mdb.models[modelName].FieldOutputRequest(
name='Ouput Request Shear_zx', 
createStepName='Isothermal linear perturbation step',
variables=('U',), 
region=regionDef,
sectionPoints=DEFAULT,
rebar=EXCLUDE)
# Deactivate Shear_zx reporting for the thermomechanical step.
mdb.models[modelName].fieldOutputRequests['Ouput Request Shear_zx'].deactivate('Thermomechanical step')	
# Delete the default Field Output Request.
del mdb.models[modelName].fieldOutputRequests['F-Output-1']

# ___________________________________
# 2. ApplyLoadsToConstraints
# ApplyLoadsToConstraintsDriver(modelName, loadVector)
loadVector = (cellVolume, cellVolume,cellVolume,cellVolume,cellVolume,cellVolume,1.0)
modelName,vecLoad = modelName, loadVector
# ___________________________________
# nodeNo = []
# for i in mdb.models[modelName].rootAssembly.allInstances[instanceName].nodes:
#     nodeNo.append(int(i.label))

# a.SetFromNodeLabels(name='CTE_part', nodeLabels=((instanceName,nodeNo),))

# region = a.sets['CTE_part']
# #                                region = a.instances[instanceName].sets['CTE_part']
# mdb.models[modelName].Temperature(name='Predefined Field-1', 
#     createStepName='Initial', region=region, distributionType=UNIFORM, 
#     crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(intemp, ))
# #                                        session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
# mdb.models[modelName].predefinedFields['Predefined Field-1'].setValuesInStep(
#                     stepName='Step-1', magnitudes=(fntemp, ))
#                     # ___________________________________

# error because cell only exist in geometry object
nodeNo = [int(i.label) for i in mdb.models[modelName].rootAssembly.allInstances[instanceName].nodes]
assembly.SetFromNodeLabels(name='CTE_part', nodeLabels=((instanceName,nodeNo),))
region = assembly.sets['CTE_part']
mdb.models[modelName].Temperature(
    name='Initial temperature 0 degree C all cells',
    createStepName='Initial',
    region=region,
    distributionType=UNIFORM,
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
    magnitudes=(0.0,))

# Create a heated-up temperature field.
mdb.models[modelName].Temperature(
    name='Temperature steady 1 degree C all cells',
    createStepName='Thermomechanical step',
    region=region,
    distributionType=UNIFORM,
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
    magnitudes=(vecLoad[6],))

# Disable the original temperature distribution in the temperature jump step.
mdb.models[modelName].predefinedFields['Initial temperature 0 degree C all cells'].resetToInitial(
    stepName='Thermomechanical step')		
    
# Create a 'Concentrated Force'	for the loadings.
Fx = vecLoad[0]
Fy = vecLoad[1]
Fz = vecLoad[2]
Shear_yz = vecLoad[3]		
Shear_zx = vecLoad[4]		
Shear_xy = vecLoad[5]		
# Setting up the shear values.
if (Shear_yz!=0.0):
    regRefPointShear_yz = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_YZ']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER SHEAR_YZ',
        createStepName='Isothermal linear perturbation step', 
        region=regRefPointShear_yz,
        cf1=Shear_yz,
        distributionType=UNIFORM,
        localCsys=None)

if (Shear_zx!=0.0):
    regRefPointShear_zx = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_ZX']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER SHEAR_ZX',
        createStepName='Isothermal linear perturbation step', 
        region=regRefPointShear_zx,
        cf1=Shear_zx,
        distributionType=UNIFORM,
        localCsys=None)

if (Shear_xy!=0.0):
    regRefPointShear_xy = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_XY']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER SHEAR_XY',
        createStepName='Isothermal linear perturbation step', 
        region=regRefPointShear_xy,
        cf1=Shear_xy,
        distributionType=UNIFORM,
        localCsys=None)

# Setting up the concentrated forces.
if (Fx!=0.0):
    regRefPointFx = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FX']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER FX', 
        createStepName='Isothermal linear perturbation step',
        region=regRefPointFx,
        cf1=Fx,
        distributionType=UNIFORM,
        localCsys=None)

if (Fy!=0.0):
    regRefPointFy = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FY']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER FY', 
        createStepName='Isothermal linear perturbation step',
        region=regRefPointFy,
        cf1=Fy,
        distributionType=UNIFORM,
        localCsys=None)

if (Fz!=0.0):
    regRefPointFz = mdb.models[modelName].rootAssembly.sets['CONSTRAINTS DRIVER FZ']
    mdb.models[modelName].ConcentratedForce(
        name='Load on CONSTRAINTS DRIVER FZ', 
        createStepName='Isothermal linear perturbation step',
        region=regRefPointFz,
        cf1=Fz,
        distributionType=UNIFORM,
        localCsys=None)

# Load cases.
if (Fx!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Pure x load',
        loads=(('Load on CONSTRAINTS DRIVER FX',1),),
        includeActiveBaseStateBC=ON)

if (Fy!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Pure y load',
        loads=(('Load on CONSTRAINTS DRIVER FY',1),),
        includeActiveBaseStateBC=ON)	

if (Fz!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Pure z load',
        loads=(('Load on CONSTRAINTS DRIVER FZ',1),),
        includeActiveBaseStateBC=ON)

if (Shear_yz!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Shear yz load',
        loads=(('Load on CONSTRAINTS DRIVER SHEAR_YZ',1),),
        includeActiveBaseStateBC=ON)

if (Shear_zx!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Shear zx load',
        loads=(('Load on CONSTRAINTS DRIVER SHEAR_ZX',1),),
        includeActiveBaseStateBC=ON)

if (Shear_xy!=0.0):
    mdb.models[modelName].steps['Isothermal linear perturbation step'].LoadCase(
        name='Shear xy load',
        loads=(('Load on CONSTRAINTS DRIVER SHEAR_XY',1),),
        includeActiveBaseStateBC=ON)

# ___________________________________




# 3. JOB
jobModelName=modelName.replace(' ','_')
mdb.Job(name=jobModelName,
    model=modelName,
    description='Applying concentrated forces at the driving points',
    type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None,
    memory=70, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
    echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
    userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, 
    multiprocessingMode=DEFAULT,
    numDomains=4, numCpus=2)


def Create_viewport(modelName):
	# Create a new viewport in which to display the model.
	unitCellViewport = session.Viewport(name=modelName)
	unitCellViewport.makeCurrent()
	unitCellViewport.maximize()


Create_viewport(modelName)
# Submit the job and open the resulting database.
mdb.jobs[jobModelName].submit(consistencyChecking=ON)
mdb.jobs[jobModelName].waitForCompletion()
fileName = jobModelName + '.odb'





























import visualization
resultODB = visualization.openOdb(path=fileName, readOnly=True)


# from CommonUCFeatures import ViewportDisplay
# ViewportDisplay(modelName,modelType, resultODB, loadVector)
modelName, rODB, lVect = modelName, resultODB, loadVector 
if (lVect!=None):
    # Ensure that each viewport links to a load case.
    # If the load case does not exist, lower the 'frame' index.
    frameIndexOffset=0
    if (lVect[0]!=0.0):
        session.Viewport(
            name='Pure X load',
            titleStyle=CUSTOM,
            customTitleString='Pure X load')
        session.viewports['Pure X load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure X load'].odbDisplay.setFrame(
            step=0,
            frame=1)
        session.viewports['Pure X load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure X load'].odbDisplay.setPrimaryVariable(
        variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S11'), )
        session.viewports['Pure X load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1
    if (lVect[1]!=0.0):
        session.Viewport(
            name='Pure Y load',
            titleStyle=CUSTOM,
            customTitleString='Pure Y load')
        session.viewports['Pure Y load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure Y load'].odbDisplay.setFrame(
            step=0,
            frame=(2-frameIndexOffset))
        session.viewports['Pure Y load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure Y load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S22'), )
        session.viewports['Pure Y load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1
    if (lVect[2]!=0.0):
        session.Viewport(
            name='Pure Z load',
            titleStyle=CUSTOM,
            customTitleString='Pure Z load')
        session.viewports['Pure Z load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure Z load'].odbDisplay.setFrame(
            step=0,
            frame=(3-frameIndexOffset))
        session.viewports['Pure Z load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure Z load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S33'), )
        session.viewports['Pure Z load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1			
    if (lVect[5]!=0.0):
        session.Viewport(
            name='Pure XY shear load',
            titleStyle=CUSTOM,
            customTitleString='Pure XY shear load')
        session.viewports['Pure XY shear load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure XY shear load'].odbDisplay.setFrame(
            step=0,
            frame=(4-frameIndexOffset))
        session.viewports['Pure XY shear load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure XY shear load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S12'), )
        session.viewports['Pure XY shear load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1			
    if (lVect[3]!=0.0):
        session.Viewport(
            name='Pure YZ shear load',
            titleStyle=CUSTOM,
            customTitleString='Pure YZ shear load')
        session.viewports['Pure YZ shear load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure YZ shear load'].odbDisplay.setFrame(
            step=0,
            frame=(5-frameIndexOffset))
        session.viewports['Pure YZ shear load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure YZ shear load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S23'), )
        session.viewports['Pure YZ shear load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1			
    if (lVect[4]!=0.0):
        session.Viewport(
            name='Pure ZX shear load',
            titleStyle=CUSTOM,
            customTitleString='Pure ZX shear load')
        session.viewports['Pure ZX shear load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Pure ZX shear load'].odbDisplay.setFrame(
            step=0,
            frame=(6-frameIndexOffset))
        session.viewports['Pure ZX shear load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Pure ZX shear load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(COMPONENT,'S13'), )
        session.viewports['Pure ZX shear load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    else:
        frameIndexOffset = frameIndexOffset + 1			
    if (lVect[6]!=0.0):
        session.Viewport(
            name='Temperature load',
            titleStyle=CUSTOM,
            customTitleString='Temperature load')
        session.viewports['Temperature load'].setValues(displayedObject=rODB)
        # Align the viewport with its load case results.
        session.viewports['Temperature load'].odbDisplay.setFrame(
            step=1,
            frame=(1))
        session.viewports['Temperature load'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        session.viewports['Temperature load'].odbDisplay.setPrimaryVariable(
            variableLabel='S',
            outputPosition=INTEGRATION_POINT,
            refinement=(INVARIANT,'Mises'), )
        session.viewports['Temperature load'].view.rotate(
            mode=MODEL,
            xAngle=270.0,
            yAngle=0.0,
            zAngle=270.0,
            drawImmediately=True)
    ####Cascade the viewports.	
    session.viewports[modelName].restore()
    session.viewports[modelName].setValues(
        origin=(10.0, -105.0),
        width=155.0,
        height=110.0)
    viewPortIndexOffset = 0
    if (lVect[0]!=0.0):
        session.viewports['Pure X load'].setValues(
            origin=(15.0, 10.0),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1
    if (lVect[1]!=0.0):
        session.viewports['Pure Y load'].setValues(
            origin=(5.0*(4-viewPortIndexOffset), 10.0-10.0*(1-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1
    if (lVect[2]!=0.0):
        session.viewports['Pure Z load'].setValues(
            origin=(5.0*(5-viewPortIndexOffset), 10.0-10.0*(2-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1
    if (lVect[5]!=0.0):
        session.viewports['Pure XY shear load'].setValues(
            origin=(5.0*(6-viewPortIndexOffset), 10.0-10.0*(3-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1
    if (lVect[3]!=0.0):
        session.viewports['Pure YZ shear load'].setValues(
            origin=(5.0*(7-viewPortIndexOffset), 10.0-10.0*(4-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1		
    if (lVect[4]!=0.0):
        session.viewports['Pure ZX shear load'].setValues(
            origin=(5.0*(8-viewPortIndexOffset), 10.0-10.0*(5-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    else:
        viewPortIndexOffset = viewPortIndexOffset + 1		
    if (lVect[6]!=0.0):
        session.viewports['Temperature load'].setValues(
            origin=(5.0*(9-viewPortIndexOffset), 10.0-10.0*(6-viewPortIndexOffset)),
            width=155.0,
            height=110.0)
    


session.viewports[modelName].setValues(
    origin=(15.0, 10.0),
    width=155.0,
    height=110.0)

volUC, lVect, rODB = cellVolume, loadVector, resultODB

scale = 1.0

# The constraints are stored in the load-vector.
# Fx, Fy, Fz, Shear_yz, Shear_zx, Shear_xy.

# Collect the various steps.
isothermalStep = rODB.steps.values()[0]
thermoMechanicalStep = rODB.steps.values()[1]
    
# The load-cases are held in frames.
# Index 0 holds the reference frame.
# 1 holds the Fx load-case.
frameFx = rODB.steps[isothermalStep.name].frames[1]
# 2 holds the Fy load-case.
frameFy = rODB.steps[isothermalStep.name].frames[2]
# 3 holds the Fx load-case.
frameFz = rODB.steps[isothermalStep.name].frames[3]
# 4 holds the Shear_xy load-case.
frameShear_xy = rODB.steps[isothermalStep.name].frames[4]
# 5 holds the Shear_yz load-case.
frameShear_yz = rODB.steps[isothermalStep.name].frames[5]
# 6 holds the Shear_zx load-case.
frameShear_zx = rODB.steps[isothermalStep.name].frames[6]

# Index 0 holds the reference frame.
frameThermomechanical = rODB.steps[thermoMechanicalStep.name].frames[1]
    
# Load case Fx.
Fx_eps0_x = frameFx.fieldOutputs['U'].values[0].data[0]
Fx_eps0_y = frameFx.fieldOutputs['U'].values[1].data[0]
Fx_eps0_z = frameFx.fieldOutputs['U'].values[2].data[0]

# Load case FY.
Fy_eps0_x = frameFy.fieldOutputs['U'].values[0].data[0]
Fy_eps0_y = frameFy.fieldOutputs['U'].values[1].data[0]
Fy_eps0_z = frameFy.fieldOutputs['U'].values[2].data[0]

# Load case Fz.
Fz_eps0_x = frameFz.fieldOutputs['U'].values[0].data[0]
Fz_eps0_y = frameFz.fieldOutputs['U'].values[1].data[0]
Fz_eps0_z = frameFz.fieldOutputs['U'].values[2].data[0]

# Load case Shear_xy.
Shear_xy_gamma0_xy = frameShear_xy.fieldOutputs['U'].values[3].data[0]

# Load case Shear_yz.
Shear_yz_gamma0_yz = frameShear_yz.fieldOutputs['U'].values[4].data[0]

# Load case Shear_zx.
Shear_zx_gamma0_zx = frameShear_zx.fieldOutputs['U'].values[5].data[0]



# Thermomechanical expansion.
Th_eps0_x = frameThermomechanical.fieldOutputs['U'].values[0].data[0]
Th_eps0_y = frameThermomechanical.fieldOutputs['U'].values[1].data[0]
Th_eps0_z = frameThermomechanical.fieldOutputs['U'].values[2].data[0]

# Material properties.
E0_x = lVect[0]/volUC/Fx_eps0_x*scale*scale
v0_xy = -Fx_eps0_y/Fx_eps0_x
v0_xz = -Fx_eps0_z/Fx_eps0_x

E0_y = lVect[1]/volUC/Fy_eps0_y*scale*scale
v0_yx = -Fy_eps0_x/Fy_eps0_y
v0_yz = -Fy_eps0_z/Fy_eps0_y

E0_z = lVect[2]/volUC/Fz_eps0_z*scale*scale
v0_zx = -Fz_eps0_x/Fz_eps0_z
v0_zy = -Fz_eps0_y/Fz_eps0_z

G0_yz = lVect[3]/volUC/Shear_yz_gamma0_yz*scale*scale

G0_zx = lVect[4]/volUC/Shear_zx_gamma0_zx*scale*scale

G0_xy = lVect[5]/volUC/Shear_xy_gamma0_xy*scale*scale

alpha0_x = Th_eps0_x/lVect[6]
alpha0_y = Th_eps0_y/lVect[6]
alpha0_z = Th_eps0_z/lVect[6]
Fx_eps0_x=Fx_eps0_x/scale/scale
Fx_eps0_y=Fx_eps0_y/scale/scale
Fx_eps0_z=Fx_eps0_z/scale/scale
Fy_eps0_x=Fy_eps0_x/scale/scale
Fy_eps0_y=Fy_eps0_y/scale/scale
Fy_eps0_z=Fy_eps0_z/scale/scale
Fz_eps0_x=Fz_eps0_x/scale/scale
Fz_eps0_y=Fz_eps0_y/scale/scale
Fz_eps0_z=Fz_eps0_z/scale/scale
Shear_yz_gamma0_yz=Shear_yz_gamma0_yz/scale/scale
Shear_zx_gamma0_zx=Shear_zx_gamma0_zx/scale/scale
Shear_xy_gamma0_xy=Shear_xy_gamma0_xy/scale/scale
	
	