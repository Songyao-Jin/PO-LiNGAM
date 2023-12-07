import numpy as np
import pandas as pd
import itertools
import networkx as nx
import warnings
import Utils_V2 as Utils
from itertools import combinations
import GIN2 


#return the combination of list
def FindCombination(Lists,N):
    return itertools.combinations(Lists,N)


#GIN test by fast HSIC
#X and Z are list, e.g., X=['X1','X2'] Z=['X3']
#Data.type=Pandas.DataFrame, where data.columns=['x1','x2',...]
def GIN(X,Z,data,alpha=0.05):
    return GIN2.GIN(X,Z,data,alpha)

    
    

def CheckOverlapBetweenAtomicUnits(data, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits, alpha):
    
    if len(activeAtomicUnits) <=1:
        return activeAtomicUnits

    #check if two atomic units have overlap.
    Sets_P = FindCombination(activeAtomicUnits, 2)
    for P in Sets_P:    
        P = list(P)
        if P[0] in newGeneratedAtomicUnits or P[1] in newGeneratedAtomicUnits:
            Z_set = []

            for unit in activeAtomicUnits:
                Z_set = Z_set + unit.surrogate_set_B
            for n in range(len(P[1].var_set)):
                Y_set = P[0].surrogate_set_A + P[1].surrogate_set_A[:n+1] 
                if GIN(Y_set,Z_set,data,alpha):    
                    numOfOverlap = len(P[1].var_set)-n
                    if len(P[0].var_set) != numOfOverlap or len(P[1].var_set) != numOfOverlap: 
                        print(P[0].var_set,'and',P[1].var_set,"have",numOfOverlap,"overlapped variables")
                        P[0].add_overlapped_units(P[1], numOfOverlap)
                        P[1].add_overlapped_units(P[0], numOfOverlap)
                        break
    

    decomposedUnits = set()
    for unit in activeAtomicUnits:
        coveredUnits = []
        for key, value in unit.overlapped_units.items():
            if len(key.var_set) == value and value <= len(unit.var_set):   
                coveredUnits.append(key)
        
        for num in range(0,len(coveredUnits)):
            Sets_P=FindCombination(coveredUnits,num+1)
            for oneSet in Sets_P:
                totalNrOfVar = 0
                disjoint = True
                for oneCoveredUnit in oneSet:
                    totalNrOfVar += unit.overlapped_units[oneCoveredUnit]
                    oneCopyedSet = list(oneSet).copy()
                    oneCopyedSet.remove(oneCoveredUnit)
                    for anotherCoveredUnit in oneCopyedSet:
                        if anotherCoveredUnit in oneCoveredUnit.overlapped_units.keys():
                            disjoint = False
                if totalNrOfVar == len(unit.var_set) and disjoint:  
                    if len(oneSet) > 1:                          
                        activeAtomicUnits = UpdateCausalRelation(oneSet,unit,activeAtomicUnits,inactiveAtomicUnits)
                        decomposedUnits.add(unit)
                        print(unit.var_set,"can be fully decomposited as",[unit.var_set for unit in oneSet])
                

    activeAtomicUnits = [x for x in activeAtomicUnits if x not in decomposedUnits]
    
    return activeAtomicUnits              
                    
         
def UpdateCausalRelation(coverUnitsSet, decompositedUnit, activeAtomicUnits, inactiveAtomicUnits):
    

    children = set()
    for unit in inactiveAtomicUnits:
        for childName in decompositedUnit.children:
            if unit.var_set == childName:
                children.add(unit)
    
    for unit in coverUnitsSet:
        unit.update_children(decompositedUnit.children)
        unit.remove_overlapped_units(decompositedUnit)
        for child in children:
            child.remove_parent(decompositedUnit.var_set) 
            child.add_parent(unit.var_set)

    return activeAtomicUnits