# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 21:48:17 2020

@author: Ellias
"""
#import regionToolset
#from abaqus import *
from abaqusConstants import *
#from caeModules import *
#from odbAccess import *
#import visualization
#import customKernel

import AbaqusUtils as au

####### INPUT DATA: Check all input data before use ########

####### Input Model Data #######
mdbName = 'C:\Work\Abaqus\TestDir/test.cae'
mdbModelName = 'TestModel'
partName = 'TestPart'
#stepName = 'Step'
borderNodeSet = 'N_Set'
#borderElementSet = 'E_Set'
rWater = 30#96.5           # water level (upperstream for upper side, downstream for down side
rFound = 5 #92             # ground level at the bottom of reservoir near dam
hw = rWater - rFound    # water depth near dam
normalArea = [1,0,0] #Norlam to water pressure face, use to define mass component cofficients (x,y,z)

SFrw = 1    # Area per node in m^2
            # Use 0 to calculate automatically according to you mesh 
            #     Note: To use automated method you should create SFM3D3 element type mesh on the surface of interest
            #     For now works correctly only with 3-node elements and X-axis normal, normalArea [1,0,0]


####### Input SP Data #######
density = 1000          # Water density
Psi = 0.73              # Water basin limit coefficient according SP

def massRzh(a):
    x = [0.0, 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7, 0.8,  0.9,  1.0 ]
    y = [0.0, 0.23, 0.36, 0.47, 0.55, 0.61, 0.66, 0.7, 0.72, 0.74, 0.74]
    return au.multiLinearApprox(x, y, a)

def Mu1Var1(a):
    x = [0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.00]
    y = [0.0,  0.22, 0.38, 0.47, 0.53, 0.57, 0.59, 0.61, 0.62, 0.63, 0.68]
    return au.multiLinearApprox(x, y, a)
    
def Mu1Var2(a):
    x = [0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.00]
    y = [0.0,  0.22, 0.35, 0.41, 0.46, 0.49, 0.52, 0.53, 0.54, 0.54, 0.55]
    return au.multiLinearApprox(x, y, a)
    
def Mu1Var3(a):
    x = [0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.00]
    y = [0.0,  0.21, 0.29, 0.35, 0.38, 0.41, 0.43, 0.44, 0.45, 0.45, 0.44]
    return au.multiLinearApprox(x, y, a)


####### END OF INPUT DATA: Check all input data before use ########


####### Main Script ##########

global myMdb
myMdb = openMdb(mdbName)

# get nodes in the node set by setname
borderNodes = myMdb.models[mdbModelName].parts[partName].sets[borderNodeSet].nodes
L = len(borderNodes)
meshNodeObj = [0] * L
myRegion = [0] * L
sumSFrw = 0
L0 = 0          #Counter for nodes with null area

#cycle thru nodes
for i in range (0,L):
    nodeCoordY = borderNodes[i].coordinates[1]
    nodeLabel = borderNodes[i].label
    
    if (nodeCoordY < rFound):
        nodeCoordY = rFound
        
    if (nodeCoordY < rWater): 
        # nodes is below water level
        if (SFrw == 0):
            SFrw = au.getNodeAreaXProj(borderNodes[i])

        if (SFrw == 0):
            L0 = L0 + 1
        else:
            #check you case by SP before use massRzh 
            mTemp=massRzh((rWater-nodeCoordY)/hw) #Надо доработать, сделав разные случаи доступными через параметр
            
            rMassPerNode = density * hw *  mTemp * Psi * SFrw 
            sumSFrw = sumSFrw + SFrw
                      
            # get node
            meshNodeObj[i] = myMdb.models[mdbModelName].parts[partName].nodes.sequenceFromLabels((nodeLabel,))
            # create region from node
            myRegion[i] = regionToolset.Region(nodes=meshNodeObj[i])
            # create engineeringFeatures Object (service)
            eFeatures = myMdb.models[mdbModelName].parts[partName].engineeringFeatures
            ## create Point Mass
            massName = 'm' + str(nodeLabel)
            eFeatures.PointMassInertia(name=massName, region = myRegion[i], mass1=normalArea[0]*rMassPerNode, mass2=normalArea[1]*rMassPerNode, mass3=normalArea[2]*rMassPerNode)
            print('Node: ' + str(nodeLabel)+' Ycoord: '+str(nodeCoordY)+' Mass: '+str(rMassPerNode))
            print ('PointMassInertia: ' + str(i+1) + ' of ' + str(L))

print ('sumSFrw = '+ str(sumSFrw))
print ('L0 = '+ str(L0))     

####### End of Main Script ########## 