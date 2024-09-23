# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024.HF1 replay file
# Internal Version: 2024_01_13-02.34.40 RELr426 190904
# Run by im_ca on Mon Sep 23 22:05:59 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=215.859359741211, 
    height=206.125)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
execfile(
    'D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_bcc_et.py', 
    __main__.__dict__)
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
#: The model "BCC_et_1" has been created.
#: Boundary conditions have been defined
#: Warning: parallelizationMethodExplicit member has been removed from Job objects.
#: Job BCC_et_1: Analysis Input File Processor completed successfully.
#: Job BCC_et_1: Abaqus/Standard completed successfully.
#: Job BCC_et_1 completed successfully. 
#: Model: D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/BCC_et_1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       2
#: Number of Node Sets:          3748
#: Number of Steps:              2
#: Fx_eps0: [406581.03, -94603.61, -94582.04]
#: Fy_eps0: [-94603.61, 406511.3, -94629.81]
#: Fz_eps0: [-94582.04, -94629.81, 406608.4]
p = mdb.models['BCC_et_1'].parts['BCC_et_1Part']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['BCC_et_1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['BCC_et_1'].parts['BCC_et_1Part']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=254.315, 
    farPlane=477.393, width=347.67, height=135.062, viewOffsetX=11.4638, 
    viewOffsetY=-7.6467)
cliCommand("""from abaqus import *""")
cliCommand("""from abaqusConstants import *""")
cliCommand("""model_name='BCC_et_1'""")
cliCommand("""eleType = DC3D8""")
cliCommand("""elemType1 = mesh.ElemType(elemCode=eleType, elemLibrary=STANDARD)""")
#* NameError: name 'mesh' is not defined
cliCommand("""import mesh""")
cliCommand("""elemType1 = mesh.ElemType(elemCode=eleType, elemLibrary=STANDARD)""")
cliCommand("""part_name = model_name+"Part\"""")
cliCommand("""model = mdb.models[model_name]""")
cliCommand("""region = model.parts[part_name].cells""")
cliCommand("""model.parts[part_name].setElementType(regions=(region,), elemTypes=(elemType1,))""")
#* ERROR: Only regions of the same dimension may be
#* selected for each element type assignment
cliCommand("""region""")
#: mdb.models['BCC_et_1'].parts['BCC_et_1Part'].cells
cliCommand("""len(region)""")
#: 0
mdb.saveAs(
    pathName='D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/bcc_et_1')
#: The model database has been saved to "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\test\bcc_et_1.cae".
session.viewports['Viewport: 1'].view.setValues(nearPlane=249.739, 
    farPlane=481.968, width=411.053, height=159.685, viewOffsetX=32.6926, 
    viewOffsetY=-16.0902)
elemType1 = mesh.ElemType(elemCode=DC3D10, elemLibrary=STANDARD)
p = mdb.models['BCC_et_1'].parts['BCC_et_1Part']
z1 = p.elements
elems1 = z1[0:11712]
pickedRegions =(elems1, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
openMdb(
    pathName='D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/bcc_et_1.cae')
#: The model database "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\test\bcc_et_1.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
a = mdb.models['BCC_et_1'].rootAssembly
a.regenerate()
a = mdb.models['BCC_et_1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
cliCommand("""import mesh""")
cliCommand("""from abaqus import *""")
cliCommand("""from abaqusConstants import *""")
cliCommand("""model_name='BCC_et_1'""")
cliCommand("""model_name='BCC_et_1'""")
p = mdb.models['BCC_et_1'].parts['BCC_et_1Part']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
cliCommand("""eleType=DC3D10""")
cliCommand("""part_name=model_name+"Part\"""")
cliCommand("""region = model.parts[part_name].elements""")
#* AccessError: mdb.models['BCC_et_1'] no longer exists
cliCommand("""part_name""")
#: 'BCC_et_1Part'
cliCommand("""region = model.parts[part_name].elements""")
#* AccessError: mdb.models['BCC_et_1'] no longer exists
cliCommand("""model_name""")
#: 'BCC_et_1'
cliCommand("""model_name""")
#: 'BCC_et_1'
cliCommand("""model_name""")
#: 'BCC_et_1'
cliCommand("""model""")
#* AccessError: mdb.models['BCC_et_1'] no longer exists
cliCommand("""model = mdb.models[model_name]""")
cliCommand("""region = model.parts[part_name].elements""")
cliCommand("""model.parts[part_name].setElementType(regions=(region,), elemTypes=(elemType1,))""")
#* 
#* The Tet shape associated with the regions
#*  does not have a valid element type.
#* 
#* Select a type that supports the Tet shape or change the
#* mesh controls to remove the presence of tet elements.
#* 
cliCommand("""region = model.parts[part_name].elements[:]""")
cliCommand("""model.parts[part_name].setElementType(regions=(region,), elemTypes=(elemType1,))""")
#* 
#* The Tet shape associated with the regions
#*  does not have a valid element type.
#* 
#* Select a type that supports the Tet shape or change the
#* mesh controls to remove the presence of tet elements.
#* 
cliCommand("""elemType1""")
#: ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
cliCommand("""elemType1""")
#: ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
cliCommand("""eleType""")
#: DC3D10
cliCommand("""elemType1 = mesh.ElemType(elemCode=eleType, elemLibrary=STANDARD)""")
cliCommand("""model.parts[part_name].setElementType(regions=(region,), elemTypes=(elemType1,))""")
