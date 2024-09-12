from abaqus import *
from abaqusConstants import *
# Import objects and functions to draw the sketch.
import sketch
# Import some GUI facilities.
from abaqus import getInputs
from abaqus import getWarningReply, YES, NO
# Used for sqrt.
import math
# # Import objects and functions allowing the definition and handling of ABAQUS parts.
import part
# Import the assembly module.
import assembly
# Import the step module.
import step
# Import the mesh module.
import mesh
# Import viewports etc.
import visualization
# Import load functions.
import load
# Import the job module.
import job
# Import the regionToolset module to handle the region tagging.
import regionToolset
# Import the database handling functions.
import odbAccess
# Import the database material handling functions.
import odbMaterial
# Import the database section handling functions.
import odbSection

import numpy as np
from typing import Dict, Any
from unit_cell import UnitCell, BCCUnitCell, FCCUnitCell, CubeUnitCell

import csv
import json


from collections import namedtuple
class AbaqusUtilities:
    '''
    Static class library 
    Provides static methods for executing abaqus scripts 
    This class contains generic methods
    '''
    @staticmethod
    def _material_definition(model_name: str, material_properties: Dict[str, Any], scale: float) :
        """
        Assume that the Mdb model has been defined,
        and only ISOTROPIC case is under consideration.

        ElasticThermalHomogenization

        Define the material properties and create material section
        """
        material_name = material_properties['name']
        mdb.models[model_name].Material(name=material_name)
        if material_properties['type'].upper() == 'ISOTROPIC':
            E = material_properties['E']
            nu = material_properties['nu']
            CTE = material_properties['CTE']
            mdb.models[model_name].materials[material_name].Elastic(
                type=ISOTROPIC,
                table=((E/(scale)**2, \
                nu),))
            mdb.models[model_name].materials[material_name].Expansion(
                type=ISOTROPIC,
                table=((CTE,),))

    @staticmethod
    def _create_section(model_name: str, section_name: str, material_name: str):
        mdb.models[model_name].HomogeneousSolidSection(
            name=section_name,
            material=material_name,
            thickness=None)

    @staticmethod
    def _assign_section(model_name: str, section_name:str):
        part_name = AbaqusUtilities._get_part_name(model_name)
        p = mdb.models[model_name].parts[part_name]
        elements = p.elements[:]
        region = p.Set(elements=elements, name='ElementSet-1')
        p.SectionAssignment(region=region, sectionName=section_name)


    @staticmethod
    def assign_material(model_name: str, material_properties: Dict[str, Any], scale: float):
        section_name = material_name = material_properties['name']
        AbaqusUtilities._material_definition(model_name, material_properties, scale)
        AbaqusUtilities._create_section(model_name, section_name, material_name)
        AbaqusUtilities._assign_section(model_name, section_name)

    @staticmethod
    def reference_points(model_name: str):
        # This function handles the creation of six reference points to act as dumShear_zx nodes 
        # for the equations set-up. Then it creates the equations.
        # With all the set defined, only the offset dimensions are needed.
        # Point holding the Fx load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Fy load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Fz load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Shear_yz load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Shear_zx load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Shear_xy load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Tx load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Ty load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))
        # Point holding the Tz load.
        mdb.models[model_name].rootAssembly.ReferencePoint(
            point=(0.0, 0.0, 0.0))

        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-1', 
            toName='Fx')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-2', 
            toName='Fy')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-3', 
            toName='Fz')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-4', 
            toName='Shear_yz')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-5', 
            toName='Shear_zx')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-6', 
            toName='Shear_xy')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-7', 
            toName='Tx')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-8', 
            toName='Ty')
        mdb.models[model_name].rootAssembly.features.changeKey(
            fromName='RP-9', 
            toName='Tz')

        keyRefPointList = mdb.models[model_name].rootAssembly.referencePoints.keys()
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[8]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER FX')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[7]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER FY')		
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[6]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER FZ')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[5]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER SHEAR_YZ')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[4]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER SHEAR_ZX')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[3]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='CONSTRAINTS DRIVER SHEAR_XY')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[2]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='Constraints Driver Tx')
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[1]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='Constraints Driver Ty')		
        refPoint = (mdb.models[model_name].rootAssembly.referencePoints[keyRefPointList[0]],)
        mdb.models[model_name].rootAssembly.Set(
            referencePoints=refPoint,
            name='Constraints Driver Tz')
    
    @classmethod
    def setup_mesh_pbc(cls, model_name: str, mesh_sens: float, dims: list):
        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        nodeset = mdb.models[model_name].rootAssembly.instances[instance_name].nodes
        nodes_coords = cls._classify_all_nodes(nodeset, mesh_sens, dims)
        nodes_matched = cls._nodes_matching(nodes_coords, mesh_sens)
        cls._create_nodeset(model_name, nodes_matched)
        cls._set_PBC(model_name, nodes_matched, dims)
        return 

    @classmethod
    def _set_PBC(cls, model_name: str, nodes_matched, dims: list):
        # The length of the whole unit cell model
        xDim, yDim, zDim = dims
 
        #  Vertex
        vertex_names = [f'Node_vertex{i}_%s'%(nodes_matched[f'vertex{i}'][0]) for i in range(1,9)] 
        vertex_equation_data = [
            {
                "name": "U2-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[1], 1),
                    (-1.0, vertex_names[0], 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)
                )
            },
            {
                "name": "U2-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[1], 2),
                    (-1.0, vertex_names[0], 2),
                    (-xDim,  'CONSTRAINTS DRIVER SHEAR_XY', 1)
                )
            },
            {
                "name": "U2-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[1], 3),
                    (-1.0, vertex_names[0], 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)
                )
            },
            {
                "name": "U3-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[2], 1),
                    (-1.0, vertex_names[0], 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)
                )
            },
            {
                "name": "U3-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[2], 2),
                    (-1.0, vertex_names[0], 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)
                )
            },
            {
                "name": "U3-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[2], 3),
                    (-1.0, vertex_names[0], 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1) 
                )
            },
            {
                "name": "U4-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[3], 1),
                    (-1.0, vertex_names[0], 1)
                )
            },
            {
                "name": "U4-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[3], 2),
                    (-1.0, vertex_names[0], 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)
                )
            },
            {
                "name": "U4-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[3], 3),
                    (-1.0, vertex_names[0], 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)
                )
            },
            {
                "name": "U5-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[4], 1),
                    (-1.0, vertex_names[0], 1)
                )
            },
            {
                "name": "U5-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[4], 2),
                    (-1.0, vertex_names[0], 2)
                )
            },
            {
                "name": "U5-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[4], 3),
                    (-1.0, vertex_names[0], 3),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)
                )
            },
            {
                "name": "U6-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[5], 1),
                    (-1.0, vertex_names[0], 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)
                )
            },
            {
                "name": "U6-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[5], 2),
                    (-1.0, vertex_names[0], 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)
                )
            },
            {
                "name": "U6-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[5], 3),
                    (-1.0, vertex_names[0], 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)
                )
            },
            {
                "name": "U7-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[6], 1),
                    (-1.0, vertex_names[0], 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)
                )
            },
            {
                "name": "U7-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[6], 2),
                    (-1.0, vertex_names[0], 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)
                )
            },
            {
                "name": "U7-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[6], 3),
                    (-1.0, vertex_names[0], 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1) 
                )
            },
            {
                "name": "U8-U1 Master node on u",
                "terms": (
                    (1.0, vertex_names[7], 1),
                    (-1.0, vertex_names[0], 1)
                )
            },
            {
                "name": "U8-U1 Master node on v",
                "terms": (
                    (1.0, vertex_names[7], 2),
                    (-1.0, vertex_names[0], 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)
                )
            },
            {
                "name": "U8-U1 Master node on w",
                "terms": (
                    (1.0, vertex_names[7], 3),
                    (-1.0, vertex_names[0], 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)
                )
            }
        ]
      
        for eq_data in vertex_equation_data:
            mdb.models[model_name].Equation(
                name=eq_data["name"],
                terms=eq_data["terms"]
            )


        # Edge
        edges = [nodes_matched[f'edge{i}'] for i in range(1, 13)]
        # edge I-IV
        for count, (i,j,k,l) in enumerate(zip(*edges[:4]), start=1):
            nodeSetLabels = [f'Node_edge{idx}_{label}' for idx, label in enumerate([i, j, k, l], start=1)]
            nodeSetLabel_I, nodeSetLabel_II, nodeSetLabel_III, nodeSetLabel_IV = nodeSetLabels
            eqEdgeName = 'U_II-U_I edge node '
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_II, 1),
                    (-1.0, nodeSetLabel_I, 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_II, 2),
                    (-1.0, nodeSetLabel_I, 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_II, 3),
                    (-1.0, nodeSetLabel_I, 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))
            # This is the U_III-U_I=FAB+FCD equation.
            eqEdgeName = 'U_III-U_I edge node '
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_III, 1),
                    (-1.0, nodeSetLabel_I, 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_III, 2),
                    (-1.0, nodeSetLabel_I, 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_III, 3),
                    (-1.0, nodeSetLabel_I, 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))
            # This is the U_IV-U_I=FCD equation.
            eqEdgeName = 'U_IV-U_I edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_IV, 1),
                    (-1.0, nodeSetLabel_I, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_IV, 2),
                    (-1.0, nodeSetLabel_I, 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_IV, 3),
                    (-1.0, nodeSetLabel_I, 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))

        # edge V-VIII
        for count, (i,j,k,l) in enumerate(zip(*edges[4:8]), start=1):
            nodeSetLabels = [f'Node_edge{idx}_{label}' for idx, label in enumerate([i, j, k, l], start=5)]
            nodeSetLabel_V, nodeSetLabel_VI, nodeSetLabel_VII, nodeSetLabel_VIII = nodeSetLabels
            # This is the U_VI-U_V=FAB equation.
            eqEdgeName = 'U_VI-U_V edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VI, 1),
                    (-1.0, nodeSetLabel_V, 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VI, 2),
                    (-1.0, nodeSetLabel_V, 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_VI, 3),
                    (-1.0, nodeSetLabel_V, 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))
            # This is the U_VII-U_V=FAB+FEF equation.
            eqEdgeName = 'U_VII-U_V edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VII, 1),
                    (-1.0, nodeSetLabel_V, 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VII, 2),
                    (-1.0, nodeSetLabel_V, 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_VII, 3),
                    (-1.0, nodeSetLabel_V, 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
            # This is the U_VIII-U_V=FEF equation.
            eqEdgeName = 'U_VIII-U_V edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VIII, 1),
                    (-1.0, nodeSetLabel_V, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_VIII, 2),
                    (-1.0, nodeSetLabel_V, 2)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_VIII, 3),
                    (-1.0, nodeSetLabel_V, 3),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

        # edge IX-XII
        for count, (i,j,k,l) in enumerate(zip(*edges[8:]), start=1):
            nodeSetLabels = [f'Node_edge{idx}_{label}' for idx, label in enumerate([i, j, k, l], start=9)]
            nodeSetLabel_IX, nodeSetLabel_X, nodeSetLabel_XI, nodeSetLabel_XII = nodeSetLabels
            # This is the U_X-U_IX=FCD equation.
            eqEdgeName = 'U_X-U_IX edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_X, 1),
                    (-1.0, nodeSetLabel_IX, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_X, 2),
                    (-1.0, nodeSetLabel_IX, 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_X, 3),
                    (-1.0, nodeSetLabel_IX, 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))
            # This is the U_XI-U_IX=FCD+FEF equation.
            eqEdgeName = 'U_XI-U_IX edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_XI, 1),
                    (-1.0, nodeSetLabel_IX, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_XI, 2),
                    (-1.0, nodeSetLabel_IX, 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_XI, 3),
                    (-1.0, nodeSetLabel_IX, 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
            # This is the U_XII-U_IX=FEF equation.
            eqEdgeName = 'U_XII-U_IX edge node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqEdgeName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_XII, 1),
                    (-1.0, nodeSetLabel_IX, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqEdgeName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_XII, 2),
                    (-1.0, nodeSetLabel_IX, 2)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqEdgeName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_XII, 3),
                    (-1.0, nodeSetLabel_IX, 3),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))

        # Face
        faces = [nodes_matched[f'face{i}'] for i in range(1, 7)]
        # face A-B
        for count, (i, j) in enumerate(zip(*faces[:2]), start=1):
            nodeSetLabel_A = f'Node_face1_{i}'
            nodeSetLabel_B = f'Node_face2_{j}'
            eqFaceName = 'U_B-U_A face node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqFaceName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_B, 1),
                    (-1.0, nodeSetLabel_A, 1),
                    (-xDim, 'CONSTRAINTS DRIVER FX', 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqFaceName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_B, 2),
                    (-1.0, nodeSetLabel_A, 2),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_XY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqFaceName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_B, 3),
                    (-1.0, nodeSetLabel_A, 3),
                    (-xDim, 'CONSTRAINTS DRIVER SHEAR_ZX', 1)))


        # This is the U_D-U_C=FCD equation.
        for count, (i, j) in enumerate(zip(*faces[2:4]), start=1):
            nodeSetLabel_C = f'Node_face3_{i}'
            nodeSetLabel_D = f'Node_face4_{j}'
            eqFaceName = 'U_D-U_C face node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqFaceName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_D, 1),
                    (-1.0, nodeSetLabel_C, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqFaceName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_D, 2),
                    (-1.0, nodeSetLabel_C, 2),
                    (-yDim, 'CONSTRAINTS DRIVER FY', 1)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqFaceName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_D, 3),
                    (-1.0, nodeSetLabel_C, 3),
                    (-yDim, 'CONSTRAINTS DRIVER SHEAR_YZ', 1)))

        # This is the U_F-U_E=FEF equation.
        for count, (i, j) in enumerate(zip(*faces[4:]), start=1):
            nodeSetLabel_E = f'Node_face5_{i}'
            nodeSetLabel_F = f'Node_face6_{j}'
            eqFaceName = 'U_F-U_E face node '
            # The 'u' equation (Abaqus coordinate 1)
            eqName = eqFaceName + str(count) + ' on u'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_F, 1),
                    (-1.0, nodeSetLabel_E, 1)))
            # The 'v' equation (Abaqus coordinate 2)
            eqName = eqFaceName + str(count) + ' on v'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0, nodeSetLabel_F, 2),
                    (-1.0, nodeSetLabel_E, 2)))
            # The 'w' equation (Abaqus coordinate 3)
            eqName = eqFaceName + str(count) + ' on w'
            mdb.models[model_name].Equation(
                name=eqName,
                terms=((1.0,nodeSetLabel_F, 3),
                    (-1.0, nodeSetLabel_E, 3),
                    (-zDim, 'CONSTRAINTS DRIVER FZ', 1)))
        # Stop the translations at Master Node 1.
        regVertex1 = mdb.models[model_name].rootAssembly.sets[vertex_names[0]]
        mdb.models[model_name].DisplacementBC(
            name='Translation stop Master Node 1', 
            createStepName='Initial',
            region=regVertex1,
            u1=SET,
            u2=SET,
            u3=SET,
            ur1=UNSET, 
            ur2=UNSET,
            ur3=UNSET,
            amplitude=UNSET,
            distributionType=UNIFORM, 
            localCsys=None)	
        print ('Boundary conditions have been defined') 

    @classmethod
    def _create_nodeset(cls, model_name: str, nodes_matched):
        '''
        create NodeSet in Abaqus.
        '''
        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        assembly = mdb.models[model_name].rootAssembly
        for name, labels in nodes_matched.items():
            for label in labels:
                assembly.SetFromNodeLabels(name=f'Node_{name}_{label}', nodeLabels=((instance_name,[label]),))
      
    @classmethod
    def _nodes_matching(cls, nodes_coords, mesh_sens):
        nodes_matched = {}
        def match_faces(face1, face2, mesh_sens, coord_indices):
            """
            Matches the nodes on both faces and returns a sorted list of labels 
            assuming there are no duplicate nodes in face1
            """
            face1_sorted = list(face1.keys())    
            face2_sorted = []

            precision = int(-np.log10(mesh_sens))

            # Build a face2 dictionary for efficient lookup
            face2_dict = {}
            for j in face2.keys():
                key = tuple(cls._round_to_precision(face2[j][index], precision) for index in coord_indices)  # 根据指定的坐标轴建立键值对
                face2_dict[key] = j

            # Matching
            for i in face1_sorted:
                key = tuple(cls._round_to_precision(face1[i][index], precision) for index in coord_indices)
                if key in face2_dict:
                    j = face2_dict[key]
                    if all(cls._is_close(face1[i][index], face2[j][index], mesh_sens) for index in coord_indices):
                        face2_sorted.append(j)

            return face1_sorted, face2_sorted

        def match_edges(edge1, edge2, mesh_sens, coord_index):
            """
            Matches nodes on two edges and returns the sorted lists of node labels.

            Parameters
            ----------
            edge1 : dict
                Dictionary of nodes on the first edge, where keys are node labels and values are coordinates.
            edge2 : dict
                Dictionary of nodes on the second edge, similar structure as edge1.
            mesh_sens : float
                Sensitivity value to determine if nodes are close enough to be considered matching.
            coord_index : int
                The index of the coordinate direction (0 for x, 1 for y, 2 for z) to compare.

            Returns
            -------
            tuple
                Two lists: one with sorted labels from edge1 and one with matching labels from edge2.
            """
            
            edge1_sorted = list(edge1.keys())
            edge2_sorted = []

            for i in edge1_sorted:
                for j in edge2.keys():
                    if cls._is_close(edge1[i][coord_index], edge2[j][coord_index], mesh_sens):  # 只比较一个指定方向
                        edge2_sorted.append(j)
                        break

            return edge1_sorted, edge2_sorted
        # Matching and sorting of faces
        face_pairs = [
        ('face1', 'face2', (1, 2)),  # (y, z) 
        ('face3', 'face4', (0, 2)),  # (x, z) 
        ('face5', 'face6', (0, 1))   # (x, y) 
        ]
        for face1_name, face2_name, coord_indices in face_pairs:
            face1, face2 = nodes_coords[face1_name], nodes_coords[face2_name]
            face1_sorted, face2_sorted = match_faces(face1, face2, mesh_sens, coord_indices)
            nodes_matched[face1_name] = face1_sorted
            nodes_matched[face2_name] = face2_sorted

        # Matching and sorting of edges (edge1-edge4: z, edge5-edge8: y, edge9-edge12: x)
        edge_groups = [
        ['edge1', ['edge2', 'edge3', 'edge4'], 2],  # (z)
        ['edge5', ['edge6', 'edge7', 'edge8'], 1],  # (y) 
        ['edge9', ['edge10', 'edge11', 'edge12'], 0]  # (x)
        ]
        for edge1_name, other_edges, coord_index in edge_groups:
            edge1 = nodes_coords[edge1_name]
            edge1_sorted = list(edge1.keys())

            for edge_name in other_edges:
                edge2 = nodes_coords[edge_name]
                _, edge2_sorted = match_edges(edge1, edge2, mesh_sens, coord_index)
                nodes_matched[edge_name] = edge2_sorted

            nodes_matched[edge1_name] = edge1_sorted  #

        # Extract vertex information
        vertex_list = ['vertex1', 'vertex2', 'vertex3', 'vertex4',
                         'vertex5', 'vertex6','vertex7','vertex8'
                        ]
        for vertex_name in vertex_list:
            nodes_matched[vertex_name] = list(nodes_coords[f'{vertex_name}'].keys())

        return nodes_matched

    @classmethod
    def _classify_all_nodes(cls, nodeset, mesh_sens: float , dims: list):
        '''
        Classifies all of the mesh nodes into one of the following categories: 
            - face node
            - edge node
            - vertex node
            - other 


        Parameters
        ----------
        nodeset : The nodeset objcect (from Abaqus)
        mesh_sens : float
            Mesh sensitivity
        dims: list
      

        Returns
        -------
        None
        '''
        xDim, yDim, zDim = (dim * 0.5 for dim in dims)
        boundaries = {
            'max_x': xDim,
            'max_y': yDim,
            'max_z': zDim,
            'min_x': -xDim,
            'min_y': -yDim,
            'min_z': -zDim,
        }

        nodes_coords = {}

        for node in nodeset:
            cls._classify_node(node, mesh_sens, boundaries, nodes_coords)

        return nodes_coords
    
    @classmethod
    def _classify_node(cls, node, mesh_sens, boundaries, nodes_coords):
        '''
        Classifies a mesh node into one of the following categories: 
            - face node
            - edge node
            - vertex node
            - other 
        The classification is based on the node's proximity to the boundaries of 
        the mesh within a certain sensitivity range.

        Parameters
        ----------
        node : The mesh node object (from Abaqus)
        boundaries : dict
            A dictionary with the following structure:
            - "min_x" : float
            - "max_x" : float
            - "min_y" : float
            - "max_y" : float
            - "min_z" : float
            - "max_z" : float
        mesh_sens : float
            Mesh sensitivity

        Returns
        -------
        None
            The function performs classification but does not return any value.

        '''
        x, y, z = node.coordinates
        label = int(node.label)

        close_to_min_x = cls._is_close(x, boundaries['min_x'], mesh_sens)
        close_to_max_x = cls._is_close(x, boundaries['max_x'], mesh_sens)
        close_to_min_y = cls._is_close(y, boundaries['min_y'], mesh_sens)
        close_to_max_y = cls._is_close(y, boundaries['max_y'], mesh_sens)
        close_to_min_z = cls._is_close(z, boundaries['min_z'], mesh_sens)
        close_to_max_z = cls._is_close(z, boundaries['max_z'], mesh_sens)

        # Classify vertices
        vertex_conditions = [
            ('vertex1', close_to_min_x and close_to_min_y and close_to_min_z),
            ('vertex2', close_to_max_x and close_to_min_y and close_to_min_z),
            ('vertex3', close_to_max_x and close_to_max_y and close_to_min_z),
            ('vertex4', close_to_min_x and close_to_max_y and close_to_min_z),
            ('vertex5', close_to_min_x and close_to_min_y and close_to_max_z),
            ('vertex6', close_to_max_x and close_to_min_y and close_to_max_z),
            ('vertex7', close_to_max_x and close_to_max_y and close_to_max_z),
            ('vertex8', close_to_min_x and close_to_max_y and close_to_max_z),
        ]

        for key, condition in vertex_conditions:
            if condition:
                if f'{key}' not in nodes_coords:
                    nodes_coords[f'{key}'] = {}
                nodes_coords[f'{key}'][label] = np.array([x, y, z])
                return

        edge_conditions = [
            ('edge1', close_to_min_x and close_to_min_y and boundaries['min_z'] < z < boundaries['max_z']),
            ('edge2', close_to_max_x and close_to_min_y and boundaries['min_z'] < z < boundaries['max_z']),
            ('edge3', close_to_max_x and close_to_max_y and boundaries['min_z'] < z < boundaries['max_z']),
            ('edge4', close_to_min_x and close_to_max_y and boundaries['min_z'] < z < boundaries['max_z']),
            ('edge5', close_to_min_x and close_to_min_z and boundaries['min_y'] < y < boundaries['max_y']),
            ('edge6', close_to_max_x and close_to_min_z and boundaries['min_y'] < y < boundaries['max_y']),
            ('edge7', close_to_max_x and close_to_max_z and boundaries['min_y'] < y < boundaries['max_y']),
            ('edge8', close_to_min_x and close_to_max_z and boundaries['min_y'] < y < boundaries['max_y']),
            ('edge9', close_to_min_y and close_to_min_z and boundaries['min_x'] < x < boundaries['max_x']),
            ('edge10', close_to_max_y and close_to_min_z and boundaries['min_x'] < x < boundaries['max_x']),
            ('edge11', close_to_max_y and close_to_max_z and boundaries['min_x'] < x < boundaries['max_x']),
            ('edge12', close_to_min_y and close_to_max_z and boundaries['min_x'] < x < boundaries['max_x']),
        ]

        for key, condition in edge_conditions:
            if condition:
                if f'{key}' not in nodes_coords:
                    nodes_coords[f'{key}'] = {}
                nodes_coords[f'{key}'][label] = np.array([x, y, z])
                return

        # Classify faces (exclude edge and vertex nodes)
        face_conditions = [
            ('face1', close_to_min_x and boundaries['min_y'] < y < boundaries['max_y'] and boundaries['min_z'] < z < boundaries['max_z']),
            ('face2', close_to_max_x and boundaries['min_y'] < y < boundaries['max_y'] and boundaries['min_z'] < z < boundaries['max_z']),
            ('face3', close_to_min_y and boundaries['min_x'] < x < boundaries['max_x'] and boundaries['min_z'] < z < boundaries['max_z']),
            ('face4', close_to_max_y and boundaries['min_x'] < x < boundaries['max_x'] and boundaries['min_z'] < z < boundaries['max_z']),
            ('face5', close_to_min_z and boundaries['min_x'] < x < boundaries['max_x'] and boundaries['min_y'] < y < boundaries['max_y']),
            ('face6', close_to_max_z and boundaries['min_x'] < x < boundaries['max_x'] and boundaries['min_y'] < y < boundaries['max_y']),
        ]

        for key, condition in face_conditions:
            if condition:
                if f'{key}' not in nodes_coords:
                    nodes_coords[f'{key}'] = {}
                nodes_coords[f'{key}'][label] = np.array([x, y, z])

    @classmethod
    def _is_close(cls, value1, value2, tolerance):
        '''
        Checks if two values are approximately equal within a specified tolerance.

        Parameters
        ----------
        value1 : float
            The first value to compare.
        value2 : float
            The second value to compare.
        tolerance : float
            Allowable difference for values to be considered close.


        Returns
        -------
        bool
            True if values are close, False otherwise.
        '''
        return abs(value1 - value2) <= tolerance

    @classmethod
    def _round_to_precision(cls, value, precision):
        '''
        Rounds a value to the specified number of decimal places.

        Parameters
        ----------
        value : float
            The value to be rounded.
        precision : int
            The number of decimal places to round to.

        Returns
        -------
        float
            The rounded value.
        '''   
        return round(value, precision)

    @staticmethod
    def _get_instance_name(model_name: str):
        return model_name + 'Instance'
    
    @staticmethod
    def _get_part_name(model_name: str):
        return model_name + 'Part'

    @classmethod
    def create_steps(cls, model_name):
        # This function allows for more complex step set-up.
        mdb.models[model_name].StaticLinearPerturbationStep(
            name='Isothermal linear perturbation step',
            description='Computation of effective elastic material properties',
            previous='Initial')
        
        mdb.models[model_name].StaticLinearPerturbationStep(
            name='Thermomechanical step',
            description='Computation of coefficients of thermal expansion',
            previous='Isothermal linear perturbation step')	
            
        # Delete the default history output request.
        # History output requests are not compatible with the load cases option.
        #mdb.models[model_name].historyOutputRequests['H-Output-1']
        mdb.models[model_name].HistoryOutputRequest(name='H-Output-1', 
        createStepName='Thermomechanical step', variables=('SJD', 'SJDA', 'SJDT', 
        'SJDTA'))
        
        # Stress only field output requests.
        mdb.models[model_name].FieldOutputRequest(
            name='Standard output request', 
            createStepName='Isothermal linear perturbation step',
            variables=('S','E','EVOL'))

        # Define an an additional field output request for the Fx RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FX']
        mdb.models[model_name].FieldOutputRequest(
            name='Output Request Fx', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        # Define an an additional field output request for the Fy RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FY']
        mdb.models[model_name].FieldOutputRequest(
            name='Ouput Request Fy', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        
        # Define an an additional field output request for the Fz RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FZ']
        mdb.models[model_name].FieldOutputRequest(
            name='Ouput Request Fz', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        # Define an an additional field output request for the Shear_xy RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_XY']
        mdb.models[model_name].FieldOutputRequest(
            name='Ouput Request Shear_xy', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        # Deactivate Shear_xy reporting for the thermomechanical step.
        
        mdb.models[model_name].fieldOutputRequests['Ouput Request Shear_xy'].deactivate('Thermomechanical step')	
        # Define an an additional field output request for the Shear_yz RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_YZ']
        mdb.models[model_name].FieldOutputRequest(
            name='Ouput Request Shear_yz', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        # Deactivate Shear_yz reporting for the thermomechanical step.
        mdb.models[model_name].fieldOutputRequests['Ouput Request Shear_yz'].deactivate('Thermomechanical step')	
        
        # Define an an additional field output request for the Shear_zx RP.
        regionDef=mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_ZX']
        mdb.models[model_name].FieldOutputRequest(
            name='Ouput Request Shear_zx', 
            createStepName='Isothermal linear perturbation step',
            variables=('U',), 
            region=regionDef,
            sectionPoints=DEFAULT,
            rebar=EXCLUDE)
        # Deactivate Shear_zx reporting for the thermomechanical step.
        mdb.models[model_name].fieldOutputRequests['Ouput Request Shear_zx'].deactivate('Thermomechanical step')	
        
        # Delete the default Field Output Request.
        del mdb.models[model_name].fieldOutputRequests['F-Output-1']

    @classmethod
    def apply_loads(cls, model_name, load_vec=(1.0,)*7):
        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        assembly = mdb.models[model_name].rootAssembly
        nodes_label = [int(i.label) for i in assembly.allInstances[instance_name].nodes]
        assembly.SetFromNodeLabels(name='CTE_part', nodeLabels=((instance_name, nodes_label),))

        region = assembly.sets['CTE_part']
        mdb.models[model_name].Temperature(
            name='Initial temperature 0 degree C all cells',
            createStepName='Initial',
            region=region,
            distributionType=UNIFORM,
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
            magnitudes=(0.0,))

        # Create a heated-up temperature field.
        mdb.models[model_name].Temperature(
            name='Temperature steady 1 degree C all cells',
            createStepName='Thermomechanical step',
            region=region,
            distributionType=UNIFORM,
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
            magnitudes=(load_vec[6],))

        # Disable the original temperature distribution in the temperature jump step.
        mdb.models[model_name].predefinedFields['Initial temperature 0 degree C all cells'].resetToInitial(
            stepName='Thermomechanical step')		
            
        # Create a 'Concentrated Force'	for the loadings.
        Fx = load_vec[0]
        Fy = load_vec[1]
        Fz = load_vec[2]
        Shear_yz = load_vec[3]
        Shear_zx = load_vec[4]
        Shear_xy = load_vec[5]
        # Setting up the shear values.
        if (Shear_yz!=0.0):
            regRefPointShear_yz = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_YZ']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER SHEAR_YZ',
                createStepName='Isothermal linear perturbation step', 
                region=regRefPointShear_yz,
                cf1=Shear_yz,
                distributionType=UNIFORM,
                localCsys=None)

        if (Shear_zx!=0.0):
            regRefPointShear_zx = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_ZX']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER SHEAR_ZX',
                createStepName='Isothermal linear perturbation step', 
                region=regRefPointShear_zx,
                cf1=Shear_zx,
                distributionType=UNIFORM,
                localCsys=None)

        if (Shear_xy!=0.0):
            regRefPointShear_xy = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER SHEAR_XY']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER SHEAR_XY',
                createStepName='Isothermal linear perturbation step', 
                region=regRefPointShear_xy,
                cf1=Shear_xy,
                distributionType=UNIFORM,
                localCsys=None)

        # Setting up the concentrated forces.
        if (Fx!=0.0):
            regRefPointFx = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FX']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER FX', 
                createStepName='Isothermal linear perturbation step',
                region=regRefPointFx,
                cf1=Fx,
                distributionType=UNIFORM,
                localCsys=None)

        if (Fy!=0.0):
            regRefPointFy = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FY']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER FY', 
                createStepName='Isothermal linear perturbation step',
                region=regRefPointFy,
                cf1=Fy,
                distributionType=UNIFORM,
                localCsys=None)

        if (Fz!=0.0):
            regRefPointFz = mdb.models[model_name].rootAssembly.sets['CONSTRAINTS DRIVER FZ']
            mdb.models[model_name].ConcentratedForce(
                name='Load on CONSTRAINTS DRIVER FZ', 
                createStepName='Isothermal linear perturbation step',
                region=regRefPointFz,
                cf1=Fz,
                distributionType=UNIFORM,
                localCsys=None)

        # Load cases.
        if (Fx!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Pure x load',
                loads=(('Load on CONSTRAINTS DRIVER FX',1),),
                includeActiveBaseStateBC=ON)

        if (Fy!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Pure y load',
                loads=(('Load on CONSTRAINTS DRIVER FY',1),),
                includeActiveBaseStateBC=ON)	

        if (Fz!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Pure z load',
                loads=(('Load on CONSTRAINTS DRIVER FZ',1),),
                includeActiveBaseStateBC=ON)

        if (Shear_yz!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Shear yz load',
                loads=(('Load on CONSTRAINTS DRIVER SHEAR_YZ',1),),
                includeActiveBaseStateBC=ON)

        if (Shear_zx!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Shear zx load',
                loads=(('Load on CONSTRAINTS DRIVER SHEAR_ZX',1),),
                includeActiveBaseStateBC=ON)

        if (Shear_xy!=0.0):
            mdb.models[model_name].steps['Isothermal linear perturbation step'].LoadCase(
                name='Shear xy load',
                loads=(('Load on CONSTRAINTS DRIVER SHEAR_XY',1),),
                includeActiveBaseStateBC=ON)

    @classmethod
    def submit_job(cls, model_name, load_vec=(1.0,)*7, display_in_viewport=True, job_params=None):
        job_name = model_name.replace(' ','_')
        mdb.Job(name=job_name,
					model=model_name,
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
        # UnitCell._create_viewport(model_name)

        mdb.jobs[job_name].submit(consistencyChecking=ON)
        mdb.jobs[job_name].waitForCompletion()

        fileName = job_name + '.odb'
        if display_in_viewport:
            resultODB = visualization.openOdb(path=fileName, readOnly=True)
            cls.viewport_display(model_name, resultODB, load_vec)
        else:
            resultODB = odbAccess.openOdb(path=fileName)
        return resultODB
# need to be delete
    @classmethod
    def get_odb(cls,model_name, load_vec=(1.0,)*7, display_in_viewport=True):
        job_name = model_name.replace(' ','_')
        fileName = job_name + '.odb'
        if display_in_viewport:
            resultODB = visualization.openOdb(path=fileName, readOnly=True)
            cls.viewport_display(model_name, resultODB, load_vec)
        else:
            resultODB = odbAccess.openOdb(path=fileName)
        return resultODB

    @staticmethod
    def _create_viewport(model_name):
        unitCellViewport = session.Viewport(name=model_name)
        unitCellViewport.makeCurrent()
        unitCellViewport.maximize()

    @classmethod
    def viewport_display(cls, model_name: str, odb_file, load_vec=None):
        '''
        odb_file:
        load_vec: loading vector

        the order in frame is "x -> y -> xy -> yz -> zx" (String)
        ''' 
        def create_viewport(name, step, frame_index, variable_label, comp1,comp2, 
        x_angle=270.0, y_angle=0.0, z_angle=270.0):
            session.Viewport(
                name=name,
                titleStyle=CUSTOM,
                customTitleString=name)
            session.viewports[name].setValues(displayedObject=odb_file)
            session.viewports[name].odbDisplay.setFrame(
                step=step,
                frame=frame_index)
            
            session.viewports[name].odbDisplay.display.setValues(
                plotState=CONTOURS_ON_DEF)
            session.viewports[name].odbDisplay.setPrimaryVariable(
                variableLabel=variable_label,
                outputPosition=INTEGRATION_POINT,
                refinement=(comp1, comp2))
            session.viewports[name].view.rotate(
                mode=MODEL,
                xAngle=x_angle,
                yAngle=y_angle,
                zAngle=z_angle,
                drawImmediately=True)

        if load_vec is None:
            return
        # name, step, variableLabel, refinement(comp1, comp2)
        ViewportData = namedtuple('ViewportData', ['name', 'step', 'variable_label', 'comp1', 'comp2'])
        viewport_data = [
            ViewportData('Pure X load', 0, 'S', COMPONENT, 'S11'),
            ViewportData('Pure Y load', 0, 'S', COMPONENT, 'S22'),
            ViewportData('Pure Z load', 0, 'S', COMPONENT, 'S33'),
            ViewportData('Shear YZ load', 0, 'S', COMPONENT, 'S23'),
            ViewportData('Shear ZX load', 0, 'S', COMPONENT, 'S13'),
            ViewportData('Shear XY load', 0, 'S', COMPONENT, 'S12'),
            ViewportData('Temperature load', 1, 'S', INVARIANT, 'Mises'),
        ]

        
        frame_index = 1
        load_case_order = [index for index, _ in sorted(enumerate(viewport_data[:-1]), key=lambda x: x[1].name)]
        for case_index in load_case_order:
            if load_vec[case_index] != 0.0:
                case_data = viewport_data[case_index]
                create_viewport(case_data.name, case_data.step, frame_index, case_data.variable_label, case_data.comp1, case_data.comp2)
                frame_index += 1
        # temperature
        if load_vec[-1] != 0.0:
            tem_data = viewport_data[-1]
            frame_index = 1
            create_viewport(tem_data.name, tem_data.step, frame_index, tem_data.variable_label, tem_data.comp1, tem_data.comp2)
       
        load_case_order = [index for index, _ in sorted(enumerate(viewport_data), key=lambda x: x[1].name)]
        # Position each viewport
        viewPortIndexOffset = 0
        for case_index in load_case_order:
            if load_vec[case_index] != 0.0:
                session.viewports[viewport_data[case_index].name].setValues(
                    origin=(15.0 + 5.0 * viewPortIndexOffset, 10.0  - 10.0 * viewPortIndexOffset),
                    width=155.0,
                    height=110.0)
                viewPortIndexOffset += 1


    @classmethod
    def calc_effective_props(cls, vol_unitcell: float, load_vec: list[float], odb_file, scale: float):
        def get_displacement(frame, indices):
            return [frame.fieldOutputs['U'].values[i].data[0] for i in indices]

        def compute_modulus(load: float, vol_unitcell: float, F_eps0: float, scale: float):
            return (load * scale * scale) / (vol_unitcell * F_eps0)

        def compute_poisson_ratio(disp_num: float, disp_denom: float):
            return -disp_num / disp_denom

        # Collect steps
        isothermalStep, thermoMechanicalStep = odb_file.steps.values()[0], odb_file.steps.values()[1]

        # Load cases are stored in frames
        frames_iso = [odb_file.steps[isothermalStep.name].frames[i] for i in range(1, 7)]
        frame_thermo = odb_file.steps[thermoMechanicalStep.name].frames[1]

        # Load case displacements
        Fx_eps0 = get_displacement(frames_iso[0], [0, 1, 2])
        Fy_eps0 = get_displacement(frames_iso[1], [0, 1, 2])
        Fz_eps0 = get_displacement(frames_iso[2], [0, 1, 2])
        print(f'Fx_eps0: {Fx_eps0}')
        print(f'Fy_eps0: {Fy_eps0}')
        print(f'Fz_eps0: {Fz_eps0}')
        # [Shear_yz, Shear_xz, Shear_xy]
        shear_gama0 = [frames_iso[i].fieldOutputs['U'].values[i].data[0] for i in range(3, 6)]

        # Thermomechanical expansion
        Th_eps0 = get_displacement(frame_thermo, [0, 1, 2])

        # Material properties
        # E0_x, E0_y, E0_z
        E0 = [compute_modulus(load_vec[i], vol_unitcell, F_eps0, scale) for i, F_eps0 in enumerate([Fx_eps0[0], Fy_eps0[1], Fz_eps0[2]])]
        
        # v0_xy, v0_xz
        # v0_yx, v0_yz
        # v0_zx, v0_zy
        v0 = [compute_poisson_ratio(Fx_eps0[1], Fx_eps0[0]), compute_poisson_ratio(Fx_eps0[2], Fx_eps0[0]),
            compute_poisson_ratio(Fy_eps0[0], Fy_eps0[1]), compute_poisson_ratio(Fy_eps0[2], Fy_eps0[1]),
            compute_poisson_ratio(Fz_eps0[0], Fz_eps0[2]), compute_poisson_ratio(Fz_eps0[1], Fz_eps0[2])]
        
        # [G0yz, G0xz, G0xy]
        G0 = [compute_modulus(load_vec[3 + i], vol_unitcell, shear_gama0[i], scale) for i in range(3)]

        alpha = [Th_eps0[i] / load_vec[6] for i in range(3)]

        # Scale displacements
        Fx_eps0 = [d / (scale * scale) for d in Fx_eps0]
        Fy_eps0 = [d / (scale * scale) for d in Fy_eps0]
        Fz_eps0 = [d / (scale * scale) for d in Fz_eps0]

        shear_gama0 = [d / (scale * scale) for d in shear_gama0]

        # Return results
        effective_props = {
            'E': E0,
            'v': v0,
            'G': G0,
            'Fx_eps0': Fx_eps0,
            'Fy_eps0': Fy_eps0,
            'Fz_eps0': Fz_eps0,
            'shear_gama0': shear_gama0,
            'alpha': alpha
        }
        
        return effective_props

      
    @staticmethod
    def save_effective_props(model_name: str, effective_props):
        E1, E2, E3 = effective_props['E']
        v12, v13, v23 = [effective_props['v'][i] for i in [0,1,3]]
        G23, G13, G12 = effective_props['G']
        alpha1, alpha2, alpha3 = effective_props['alpha']

        props_data ={
            'E1': E1, 'E2': E2, 'E3': E3,
            'v12': v12, 'v13': v13, 'v23': v23,
            'G23': G23, 'G13': G13, 'G12': G12,
            'alpha1': alpha1, 'alpha2': alpha2, 'alpha3': alpha3
        }

        with open(f'{model_name}_effective_props.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(props_data.keys()) 
            writer.writerow(props_data.values())

    @classmethod
    def create_report_json(cls, model_name: str, unit_cell: UnitCell, effective_props, scale: float):
        file_name = model_name + '.rpt'

        converted_props = AbaqusUtilities._convert_float32_to_float64(effective_props)


        data = {
            'scale': scale,
            'mesh_params': unit_cell.get_mesh_params(),
            'geo_params_scaled': unit_cell.get_geo_params(),
            'material_properties': unit_cell.get_material_properties(),
            'effective_props': converted_props
        }

        with open(file_name, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4, ensure_ascii=False)

    @staticmethod
    def _convert_float32_to_float64(data):
        if isinstance(data, np.ndarray): 
            if data.dtype == np.float32: 
                return data.astype(np.float64)  
        elif isinstance(data, list): 
            return [AbaqusUtilities._convert_float32_to_float64(item) if isinstance(item, np.float32) else item for item in data]
        elif isinstance(data, dict):
            return {key: AbaqusUtilities._convert_float32_to_float64(value) for key, value in data.items()}
 
        elif isinstance(data, np.float32): 
            return np.float64(data) 
        return data


class BCCAbaqusUtilities(AbaqusUtilities):
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
        radiusInter = geo_params_scaled['radiusInter']
        radiusExter_x, radiusExter_y, radiusExter_z = (geo_params_scaled[k] for k in ['radiusExter_x', 'radiusExter_y', 'radiusExter_z'])

        Mdb()
        model = mdb.Model(name=model_name)
        if(model_name != 'Model-1'):
            del mdb.models['Model-1']

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

        target_vector = np.array([xDimension, yDimension, zDimension])
        initial_vector = np.array([0, 0, zDimension]) 
        rotation_axis = np.cross(initial_vector, target_vector)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        rotation_angle = degrees(acos(np.dot(initial_vector, target_vector) / 
                            (np.linalg.norm(initial_vector) 
                            * np.linalg.norm(target_vector))))
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


        # ______________________________ #
        #      meshing 
        # ______________________________ #
        deviationFactor = mesh_params_scaled['deviationFactor']
        minSizeFactor = mesh_params_scaled['minSizeFactor']
        meshSize = mesh_params_scaled['meshSize']
      
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


        part_name = AbaqusUtilities._get_part_name(model_name)
        partMeshBCC = model.PartFromMeshMirror(name=part_name,part=partMeshBCC_1_2,
        point1=(0,yDimension-radiusExter_x,zDimension-radiusExter_x),
        point2=(-1,yDimension-radiusExter_x,zDimension-radiusExter_x))

        instance_name = AbaqusUtilities._get_instance_name(model_name) 
        BCCInst = assembly.Instance(name=instance_name, part=partMeshBCC, dependent=ON)

        del model.parts['BCC_Cell_1_8']
        del model.parts['MeshPartBCC_1_8']
        del model.parts['MeshPartBCC_1_4']
        del model.parts['MeshPartBCC_1_2']

class CubeAbaqusUtilities(AbaqusUtilities):
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

        elemType1 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
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

class FCCAbaqusUtilities(AbaqusUtilities):
    pass

class OCTAbaqusUtilities(AbaqusUtilities):
    pass
        
class AbaqusUtilityFactory:
    @staticmethod
    def get_utilities(unit_cell: UnitCell):
        if isinstance(unit_cell, BCCUnitCell):
            return BCCAbaqusUtilities()
        if isinstance(unit_cell, CubeUnitCell):
            return CubeAbaqusUtilities()
        elif isinstance(unit_cell, FCCUnitCell):
            return FCCAbaqusUtilities()
        else:
            raise AbaqusUtilities()



        


        





        

        


    
       



    