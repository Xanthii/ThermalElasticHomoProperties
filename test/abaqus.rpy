# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024.HF1 replay file
# Internal Version: 2024_01_13-02.34.40 RELr426 190904
# Run by im_ca on Fri Sep 13 16:17:36 2024
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
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
execfile(
    'D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py', 
    __main__.__dict__)
#* ImportError: No module named 'unit_cell'
#* File 
#* "D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py", 
#* line 7, in <module>
#*     from lattice_homo.problem import ElasticThermalStrategy
#* File 
#* "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\src\lattice_homo\problem.py", 
#* line 2, in <module>
#*     from unit_cell import UnitCell
execfile(
    'D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py', 
    __main__.__dict__)
#* ImportError: No module named 'unit_cell'
#* File 
#* "D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py", 
#* line 7, in <module>
#*     from lattice_homo.problem import ElasticThermalStrategy
#* File 
#* "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\src\lattice_homo\problem.py", 
#* line 2, in <module>
#*     from unit_cell import UnitCell
execfile(
    'D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py', 
    __main__.__dict__)
#* ImportError: No module named 'unit_cell'
#* File 
#* "D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py", 
#* line 7, in <module>
#*     from lattice_homo.problem import ElasticThermalStrategy
#* File 
#* "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\src\lattice_homo\problem.py", 
#* line 3, in <module>
#*     from lattice_homo.abaqus_utilities import AbaqusUtilityFactory
#* File 
#* "D:\OneDrive\Code\Toplogy\Homogenization\ThermalElasticHomoProperties\src\lattice_homo\abaqus_utilities.py", 
#* line 35, in <module>
#*     from unit_cell import UnitCell, BCCUnitCell, FCCUnitCell, CubeUnitCell
execfile(
    'D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/test_cube.py', 
    __main__.__dict__)
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
#: The model "Cube_comparision" has been created.
#: Boundary conditions have been defined
#: Warning: parallelizationMethodExplicit member has been removed from Job objects.
#: Job Cube_comparision: Analysis Input File Processor completed successfully.
#: Job Cube_comparision: Abaqus/Standard completed successfully.
#: Job Cube_comparision completed successfully. 
#: Model: D:/OneDrive/Code/Toplogy/Homogenization/ThermalElasticHomoProperties/test/Cube_comparision.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       2
#: Number of Node Sets:          874
#: Number of Steps:              2
#: Fx_eps0: [100.0, -30.0, -30.0]
#: Fy_eps0: [-30.0, 100.0, -30.0]
#: Fz_eps0: [-30.0, -30.0, 100.0]
p = mdb.models['Cube_comparision'].parts['Cube_comparisionPart']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Temperature load'].setValues(origin=(316.171844482422, 
    12.75))
session.viewports['Shear ZX load'].setValues(origin=(211.406234741211, 
    4.72222137451172))
session.viewports['Shear YZ load'].setValues(origin=(114.843742370605, 
    -9.20833587646484))
session.viewports['Shear XY load'].setValues(height=106.25)
session.viewports['Shear XY load'].setValues(origin=(92.1093673706055, 
    -7.31945037841797))
session.viewports['Pure Z load'].setValues(origin=(44.5312461853027, 
    -6.84722137451172))
session.viewports['Pure Y load'].setValues(origin=(53.6718711853027, 
    2.83333587646484))
session.viewports['Pure X load'].setValues(origin=(37.9687461853027, 
    17.4722213745117))
