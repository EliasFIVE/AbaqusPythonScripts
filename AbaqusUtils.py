# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:22:09 2020

@author: EliasFive
"""
import regionToolset
from abaqusConstants import *
import math


# Use for import testing
def importTestAbaqusUtils (printedString="Q"):
        print(" Import is OK:_" + printedString)
        
def multiLinearApprox(x, y, a):
    L = len(x)
    i = 0
    if (L != len(y)):
        print ('X and Y arrays are not correlated')
        return -1
    if (a < x[0] or a > x[L-1]):
        print ('parameter out of defined range')
        return -1
    while (i < L):
        if (a>=x[i] and a<=x[i+1]):
            return linearApprox(x[i], x[i+1], y[i], y[i+1], a)
        i = i+1
    print ('Error: value not found')
    return -1
    
def linearApprox(x0, x1, y0, y1, a):
    return (a-x0) / (x1-x0) * (y1 - y0) + y0


def getNodeArea(node):
    nodeElements3n = []
    nodeElements = node.getElements()
    elementArea = 0
    elementCorners = 3
    for nodeElement in nodeElements:
        if (str(nodeElement.type) == 'SFM3D3' or str(nodeElement.type) == 'S3R'):
            nodeElements3n.append(nodeElement)
    for nodeElement3n in nodeElements3n:
        elementArea = elementArea + getElementArea(nodeElement3n)
    if (len(nodeElements3n) != 0):
        nodeArea          = elementArea         / elementCorners
    else:
        nodeArea = 0
    print (nodeArea)
    return nodeArea
    
def getNodeAreaXProj(node):
    nodeElements3n = []
    nodeElements = node.getElements()
    elementArea = 0
    elementCorners = 3
    for nodeElement in nodeElements:
        if (str(nodeElement.type) == 'SFM3D3' or str(nodeElement.type) == 'S3R'):
            nodeElements3n.append(nodeElement)
    for nodeElement3n in nodeElements3n:
        elementArea = elementArea + getElementAreaXProj(nodeElement3n)
    if (len(nodeElements3n) != 0):
        nodeArea          = elementArea         / elementCorners
    else:
        nodeArea = 0
    print (nodeArea)
    return nodeArea

def getElementArea(element):
    # element area = 1/2 * A(x)B, 
    nodes = element.getNodes()
    a0= nodes[1].coordinates[0] - nodes[0].coordinates[0]
    a1 = nodes[1].coordinates[1] - nodes[0].coordinates[1]
    a2 = nodes[1].coordinates[2] - nodes[0].coordinates[2]
    b0 = nodes[2].coordinates[0] - nodes[0].coordinates[0]
    b1 = nodes[2].coordinates[1] - nodes[0].coordinates[1]
    b2 = nodes[2].coordinates[2] - nodes[0].coordinates[2]
    x = a1*b2-a2*b1
    y = a2*b0-a0*b2
    z = a0*b1-a1*b0
    S = 0.5 * math.sqrt(x*x + y*y + z*z)
    if (len(nodes) == 4): 
        # 4-node shell or surface element
        S = S * 2.0
    return S

def getElementAreaXProj(element):
    nodes = element.getNodes()
    a1 = nodes[1].coordinates[1] - nodes[0].coordinates[1]
    a2 = nodes[1].coordinates[2] - nodes[0].coordinates[2]
    b1 = nodes[2].coordinates[1] - nodes[0].coordinates[1]
    b2 = nodes[2].coordinates[2] - nodes[0].coordinates[2]
    x = a1*b2-a2*b1
    S = 0.5 * math.sqrt(x * x)
    if (len(nodes) == 4): 
        S = S * 2.0
    return S