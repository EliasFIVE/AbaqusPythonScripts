import regionToolset
from abaqusConstants import *
from caeModules import *
from odbAccess import *
import numpy as np

#Double asymptotic bound for seismic analyses

#To use this approach, you should: 
#   1) Create  a moded with an extended base area 
#       and perform a calculation for a single impact 
#   2) copy the main model and configure the transfer of movements 
#       at the border using the submodeling approach
#   3) perform a calculation for submodel
#   4) configure the main model for analysis
#   5) !!!! SAVE *.cae !!!! this is very important
#   6) Run this script
#   7) Check created eFeatures (Springs/Dashpots) in assembly of you main model
#   8) Save *.cae just in case
#   9) Now use this model however you want

#Only uniform foundatian version
#Only regular mesh version
#Only 2D model version for now

#=====================User Input======================================
#Data names
mdbName = 'Dam2D.cae'
mdbModelName = 'Foundation10000M' #name of target model for creation of BC
partName = 'Foundation10000' 
borderNodeSet = 'BorderNodes'
borderOrientation = 'vertical' #vertical, horizontal, sphere 
odbName = 'SubFoundation10000.odb' #Calcultation result of submodel with unit load and BC from extended faundation
stepName = 'Static'

#multidirectional impacts on extendet model
#you shoul apply sepparate calculations for ol directoins
multidirectional = True
stepNameX = stepName + 'X'
stepNameY = stepName + 'Y'

#MaterialProps
uniformFaundation = True
poissonRatio = 0.3
elasticModulus = 26000000000
density = 2500

#MeshProps
uniformMesh = True
nodeDistance = 400
#=====================End of user Input======================================
 
#=====================Utils======================================
def GetSubsetFieldValues (fieldOutputName,frame,nodeSetName,csys):
    field = frame.fieldOutputs[fieldOutputName]
    setField = field.getSubset(region=nodeSetName)
    if (csys != 0):
        transformedSetField = setField.getTransformedField(datumCsys=csys)
    else:
        transformedSetField = setField
    return transformedSetField.values

def ValueCombination(Xval,Yval):
    valueX =  Xval.data
    valueY =  Yval.data
    value = [valueX[0],valueY[1]]
    return value

#=====================Main=======================================

# In case of nonuniform foundation damping patamets should calculated for each node separately according to material props
if (uniformFaundation):
    shearModulus  = elasticModulus/(2*(1+poissonRatio))
    Vs = pow(shearModulus/density,0.5)
    s = 1 - 2*poissonRatio/(2*(1 - poissonRatio))
    Vc = Vs/s
    dampN = density*Vc
    dampT = density*Vs
    if(borderOrientation == 'vertical'):
        damp = np.array([dampN, dampT])     #damping for vertical plane border
    elif(borderOrientation == 'horizontal'):
        damp = np.array([dampT, dampN])     #damping for horizontal plane border
    elif(borderOrientation == 'sphere'):
        print('Sphere border not applyed, set damp=0')  #damping for sphere border
        damp = 0
    else:
        print('ERROR! border type not set, or set incorrectly')

# In case of nonuniform mesh this patamets should calculated for each node separately
if(uniformMesh):
    if(nodeDistance == 0):
        print('ERROR! zero node distance set for uniform mesh')
    dampNodal = damp*nodeDistance

# Read database
myMdb = openMdb(mdbName)
myOdb = visualization.openOdb(path=odbName)

# set all nodes in a instance to variable
instanceNodes = myMdb.models[mdbModelName].rootAssembly.instances[partName + '-1'].nodes
# get nodes in the node set by setname
meshNodes = myMdb.models[mdbModelName].parts[partName].sets[borderNodeSet].nodes
odbMeshInst = myOdb.rootAssembly.instances[partName.upper()+'-1'].nodeSets[borderNodeSet.upper()]
odbMeshNodes = odbMeshInst.nodes
print('meshNodes len : ' + str(len(meshNodes)))
print('odbMeshNodes len : ' + str(len(odbMeshNodes)))

#read odb data for U and RF
if (multidirectional):
    myStepX = myOdb.steps[stepNameX] 
    myStepY = myOdb.steps[stepNameY]
    valuesUX = GetSubsetFieldValues ('U',myStepX.frames[-1],odbMeshInst,0)
    valuesUY = GetSubsetFieldValues ('U',myStepY.frames[-1],odbMeshInst,0)
    valuesRFX = GetSubsetFieldValues ('RF',myStepX.frames[-1],odbMeshInst,0)
    valuesRFY = GetSubsetFieldValues ('RF',myStepY.frames[-1],odbMeshInst,0)

else:   
    myStep = myOdb.steps[stepName] 
    valuesU = GetSubsetFieldValues ('U',myStep.frames[-1],odbMeshInst,0)
    valuesRF = GetSubsetFieldValues ('RF',myStep.frames[-1],odbMeshInst,0)

# prepare Dictionarys
dictionaryUValuesByLabel = {}
if (multidirectional):
    for i in range (len(valuesUX)):
        label = valuesUX[i].nodeLabel
        data = ValueCombination(valuesUX[i],valuesUY[i])
        dictionaryUValuesByLabel[label] = data
else:
    for i in range (len(valuesU)):
        label = valuesU[i].nodeLabel
        data = valuesU[i].data
        dictionaryUValuesByLabel[label] = data 
print('dictionaryUValuesByLabel created with len: ' + str(len(dictionaryUValuesByLabel)))

dictionaryRFValuesByLabel = {}
if (multidirectional):
    for i in range (len(valuesRFX)):
        label = valuesRFX[i].nodeLabel
        data = ValueCombination(valuesRFX[i],valuesRFY[i])
        dictionaryRFValuesByLabel[label] = np.array(data)
else:
    for i in range (len(valuesRF)):
        label = valuesRF[i].nodeLabel
        data = valuesRF[i].data
        dictionaryRFValuesByLabel[label] = np.array(data) 
print('dictionaryRFValuesByLabel created with len: ' + str(len(dictionaryRFValuesByLabel)))

#not used
# nodes from model 
dictionaryMeshNodesByLabel = {
    node.label: node
    for node in meshNodes
}
print('dictionaryBorderMeshNodesByLabel created with len: ' + str(len(dictionaryMeshNodesByLabel)))

#not used
#nodes from odb
dictionaryOdbMeshNodesByLabel = {    
    node.label: node
    for node in odbMeshNodes
}
print('dictionaryBorderMeshNodesByLabel created with len: ' + str(len(dictionaryOdbMeshNodesByLabel)))

#not used
#Node sets match control
dictionaryMeshNodesByLabelSet = set(dictionaryMeshNodesByLabel)
dictionaryOdbMeshNodesByLabelSet = set(dictionaryOdbMeshNodesByLabel)
print('Node sets match control:')
print(dictionaryMeshNodesByLabelSet==dictionaryOdbMeshNodesByLabelSet)

#calculate spring stiffness
dictionaryStiffnesByLabel = {}
for node in meshNodes:
    label = node.label
    value = abs(dictionaryRFValuesByLabel[node.label]/dictionaryUValuesByLabel[node.label])
    dictionaryStiffnesByLabel[label] = value

print('dictionaryStiffnesByLabel created with len: ' + str(len(dictionaryStiffnesByLabel)))

#stiffnes control
def StiffnesTest2D(labelInt):
    print(dictionaryRFValuesByLabel[labelInt])
    print(dictionaryUValuesByLabel[labelInt])
    print(dictionaryStiffnesByLabel[labelInt])
    a1 = abs(dictionaryRFValuesByLabel[labelInt][0]/dictionaryUValuesByLabel[labelInt][0]) - dictionaryStiffnesByLabel[labelInt][0]
    a2 = abs(dictionaryRFValuesByLabel[labelInt][1]/dictionaryUValuesByLabel[labelInt][1]) - dictionaryStiffnesByLabel[labelInt][1]
    if(a1 == 0):
        print('Component 1: Correct')
    else:
        print ('Component 1: error ' + str(a1))
    if(a2 == 0):
        print ('Component 2: Correct')
    else:
        print('Component 2: error '+ str(a2))
    return

#print('Stiffnes control:')
#StiffnesTest2D(1)

dictionaryDampByLabel = {}
#calculate damping by nodal material
#scal damping by nodal area
#fill in the dictionary
for node in meshNodes:
    dictionaryDampByLabel[node.label] = dampNodal
print('dictionaryDampByLabel created with len: ' + str(len(dictionaryDampByLabel)))

#create eFeatures
eFeatureList = []
for node in meshNodes:
    meshNodeObj = instanceNodes.sequenceFromLabels((node.label,))
    region = regionToolset.Region(nodes=meshNodeObj)
    
    # get stiffnes
    KX = dictionaryStiffnesByLabel[node.label][0]
    KY = dictionaryStiffnesByLabel[node.label][1]
    
    # create engineeringFeatures Object (service)
    sFeatures = myMdb.models[mdbModelName].rootAssembly.engineeringFeatures
    
    # create spring-dashpot
    springName1 = 'Spr1N' + str(node.label)
    springName2 = 'Spr2N' + str(node.label)    
    spring1 = sFeatures.SpringDashpotToGround(name=springName, region = region, dof=1, orientation=None, springBehavior=ON, dashpotBehavior=OFF, springStiffness=KX)
    spring2 = sFeatures.SpringDashpotToGround(name=springName, region = region, dof=2, orientation=None, springBehavior=ON, dashpotBehavior=OFF, springStiffness=KY)
    
    eFeatureList.append(spring2)
    eFeatureList.append(spring1)

#create dashpots
for node in meshNodes:
    meshNodeObj = instanceNodes.sequenceFromLabels((node.label,))
    region = regionToolset.Region(nodes=meshNodeObj)

    # get stiffnes
    DX = dictionaryDampByLabel[node.label][0]
    DY = dictionaryDampByLabel[node.label][1]

    # create engineeringFeatures Object (service)
    dFeatures = myMdb.models[mdbModelName].rootAssembly.engineeringFeatures
    
    # create spring-dashpot
    dashpotName1 = 'Dsh1N' + str(node.label)
    dashpotName2 = 'Dsh2N' + str(node.label)
    if(borderOrientation == 'sphere'): 
        print('Sphere border not applyed')
    else:   
        dashpot1 = dFeatures.SpringDashpotToGround(name=dashpotName, region = region, dof=1, orientation=None, springBehavior=OFF, dashpotBehavior=ON, dashpotCoefficient=DX)
        dashpot2 = dFeatures.SpringDashpotToGround(name=dashpotName, region = region, dof=2, orientation=None, springBehavior=OFF, dashpotBehavior=ON, dashpotCoefficient=DY)
    
    eFeatureList.append(dashpot1)    
    eFeatureList.append(dashpot2)

#def ClearSpringDashpots():
#    for eFeature in eFeatureList:
#    eFeatureList.clear()    
    
    