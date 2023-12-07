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

    


def IdentifyCausalStructureAmongKnownAtomicUnits(data, activeAtomicUnits, inactiveAtomicUnits, alpha, maxNrOfParentUnits):
    index = list(activeAtomicUnits)  
    
    nrOfFound = 1   
    while nrOfFound != 0 and len(activeAtomicUnits) > 1:
        nrOfFound = 0   
    
        activeAtomicUnits_copy = activeAtomicUnits.copy()
        for unit in activeAtomicUnits_copy: 
            Z_set = activeAtomicUnits.copy()
            Z_set.remove(unit)
            
            for key in unit.overlapped_units.keys():
                if key in Z_set:
                    Z_set.remove(key)
            
            unit_surrogate_B = unit.surrogate_set_B
            Z_set_surrogate_A = list()
            for z in Z_set:
                Z_set_surrogate_A += z.surrogate_set_A

            foundParentsFlag = False
            for NrPa in range(1,maxNrOfParentUnits+1):
                Sets_P = FindCombination(Z_set,NrPa)
                for P in Sets_P:    
                    Y_set_surrogate_B = list()
                    for y in P:
                        Y_set_surrogate_B += y.surrogate_set_B   
                    
                    if GIN(Y_set_surrogate_B, Z_set_surrogate_A, data, alpha) and NrPa>1: 
                        size = len(Y_set_surrogate_B)
                        mark =False
                        for num in range(size):
                            Sets_Q = FindCombination(Y_set_surrogate_B,num+1)
                            for Q in Sets_Q:
                                new_Y_set_surrogate_B = [x for x in Y_set_surrogate_B if x not in Q] 
                                if not GIN(new_Y_set_surrogate_B, Z_set_surrogate_A, data, alpha): 
                                    Y_set_surrogate_B = new_Y_set_surrogate_B
                                    mark = True
                                    break
                            if mark == True:
                                break
                    Y_set_surrogate_B.append(unit_surrogate_B[0])       

                    
                    if GIN(Y_set_surrogate_B, Z_set_surrogate_A, data, alpha): 
                        print(Y_set_surrogate_B, Z_set_surrogate_A, "satisfies GIN")
                        for y in P:
                            unit.add_parent(y.var_set)
                            y.add_child(unit.var_set)
                        print(unit.var_set,"has parents:",unit.parents)
                        inactiveAtomicUnits.append(unit)
                        activeAtomicUnits.remove(unit)
                        
                        for changedParentUnit in P:
                            activeAtomicUnits,inactiveAtomicUnits = TryFindSmallerUnits(changedParentUnit, activeAtomicUnits,inactiveAtomicUnits, None, data, alpha)
                        
                        foundParentsFlag = True
                        nrOfFound = nrOfFound+1
                        break
                if foundParentsFlag:
                    break
            if nrOfFound != 0:
                break
    
    
    return  activeAtomicUnits, inactiveAtomicUnits    
                    
                
    
        

#Since changedParentUnit has updated children, check if changedParentUnit includes enough good children in activeatomicUnit and its overlapped units to subdivide them into sub units.
def TryFindSmallerUnits(changedParentUnit, activeAtomicUnits,inactiveAtomicUnits, currentNumOfLatentVariables, data, alpha):
    
    containSmallUnit = False
    
    allActiveOverlapUnits = set()
    for overlapAct in changedParentUnit.overlapped_units.keys():
        if overlapAct in activeAtomicUnits:
            allActiveOverlapUnits.add(overlapAct)
    allActiveOverlapUnits.add(changedParentUnit)
    allOverlapUnitsName=set()
    for unit in allActiveOverlapUnits:
        allOverlapUnitsName.add(tuple(unit.var_set))
    

    childCollection = set()  
    otherChildCollection = set()
    for inactUnit in inactiveAtomicUnits:
        for overlapUnit in allActiveOverlapUnits:
            if inactUnit.var_set in overlapUnit.children:
                parents_set = set(tuple(parent) for parent in inactUnit.parents)
                if parents_set.issubset(allOverlapUnitsName):
                    childCollection.add(inactUnit)  
                else:
                    otherChildCollection.add(inactUnit)      
    
    

    Grlen = 2
    clusterList = []   
    maxAtomicUnitSize = len(changedParentUnit.var_set)-1
    while(Grlen-1 <= maxAtomicUnitSize):
        surrogate_childCollection = []
        for unit in childCollection:
            surrogate_childCollection.extend(unit.surrogate_set_A)


        
        Sets_P=FindCombination(surrogate_childCollection,Grlen) 
        for Y_set_surrogate in Sets_P:
            Y_set_surrogate = list(Y_set_surrogate)
            units=set()    

            for act in childCollection:
                for var in Y_set_surrogate:
                    if var in act.surrogate_set_A:
                        units.add(act) 
                            
            if len(units)>1 and units not in clusterList:
                Z_set = childCollection.copy()
                Z_set_surrogate = list()
                # print([unit.var_set for unit in units])
                for at in list(units):    
                    try:
                        Z_set.remove(at)
                    except ValueError:
                        print("The unit",at.var_set,"may missed from Z_set.")
                        print("The original Z_set is,",[unit.var_set for unit in childCollection])
                        print("Current Z_set is",[unit.var_set for unit in Z_set])
                        print("The list of units that need to be removed is",[unit.var_set for unit in units])
                    
                    for key, value in at.overlapped_units.items():
                        if len(key.var_set) == value and value <= len(at.var_set):   
                            if key in Z_set:
                                Z_set.remove(key)
                    
                for z in Z_set:
                    Z_set_surrogate.append(z.surrogate_set_A[0])
                small_Y_set_surrogate = Y_set_surrogate.copy()
                small_Y_set_surrogate.pop()
                if GIN(Y_set_surrogate,Z_set_surrogate,data,alpha) and not GIN(small_Y_set_surrogate,Z_set_surrogate,data,alpha)  :    
                    clusterList.append(units)   
                    print(Y_set_surrogate,Z_set_surrogate, "satisfies GIN")
    
    
        # Find pure element and check if they can form a new atomic variable set
        mergedClusterList = []
        for baseCluster in clusterList:
            if baseCluster not in mergedClusterList: 
                for unit in baseCluster:
                    pureChildsetA = set() 
                    base = baseCluster.copy()
                    base.remove(unit)     
                          
                    surrogate_pureChildSetB =[] 
                    for element in base:
                        surrogate_pureChildSetB.extend(element.surrogate_set_A)
                    if len(surrogate_pureChildSetB) >= Grlen-1:   
                        for cluster in clusterList:
                            if cluster not in mergedClusterList:
                                difference_unit = cluster.difference(base)
       
                                if len(cluster)==len(baseCluster) and len(difference_unit) == 1:   
                                    pureChildsetA.update(difference_unit)

                        
                        if len(pureChildsetA) >= (Grlen-1):

                            good_base = True
                            for unit in base:
                                set_as_list = list(pureChildsetA)
                                clu = set([unit])
                                clu.update(set_as_list[0:Grlen - len(unit.var_set)]) 
                                if clu not in clusterList:
                                    good_base = False   
                                    
                            if good_base:  
                                fullPureChildren = pureChildsetA.union(base)  
                                for cluster in clusterList:
                                    if cluster.issubset(fullPureChildren):
                                        mergedClusterList.append(cluster) 
                                print("allPureChildren: ", [x.var_set for x in fullPureChildren])
                                print("mergedClusterList: ",len(mergedClusterList))
                                
                                
                                ignore = False
                                for unit in allActiveOverlapUnits:
                                    unit_set = set(tuple(child) for child in unit.children)
                                    fullPureChildren_Name = set(tuple(child.var_set) for child in fullPureChildren)
                                    if (Grlen-1) == len(unit.var_set) and fullPureChildren_Name.issubset(unit_set):
                                        ignore = True
                                if not ignore:
                                    containSmallUnit = True
                                    print("")
                                    print("Need to recalculate the unit set:", [unit.var_set for unit in allActiveOverlapUnits],", because the children:",[unit.var_set for unit in fullPureChildren])
                                    break
                
                    if containSmallUnit:          
                        break              
                if containSmallUnit:
                    break              
        
        if containSmallUnit:
            break   
        Grlen += 1    
        
    
    allUnits = activeAtomicUnits + inactiveAtomicUnits

    if containSmallUnit: 
        for otherC in otherChildCollection:
            trickyChild = False
            for inact in inactiveAtomicUnits:
                if inact not in otherChildCollection and inact not in childCollection and otherC.var_set in inact.children:
                    trickyChild = True
            if not trickyChild:    
                for unit in allUnits:
                    if otherC.var_set in unit.children:
                        unit.remove_child(otherC.var_set)
                if otherC in inactiveAtomicUnits:
                    inactiveAtomicUnits.remove(otherC)
                activeAtomicUnits.append(otherC)  
            else:
                otherC.parent=[]
                for al in allUnits:
                    if otherC.var_set in al.parents:
                        al.remove_child(otherC.var_set)                  
        
        
        for child in childCollection:
            for al in allActiveOverlapUnits:
                if al.var_set in child.parents:
                    child.remove_parent(al.var_set)
            if child in inactiveAtomicUnits:
                inactiveAtomicUnits.remove(child)
        activeAtomicUnits = activeAtomicUnits + list(childCollection)            
        
        for unit in allActiveOverlapUnits:
            if unit in activeAtomicUnits and "x" not in unit.var_set[0]:
                activeAtomicUnits.remove(unit)
            for al in allUnits:
                if unit in al.overlapped_units.keys():
                    al.overlapped_units.pop(unit)
               
        print("Current Active Atomic Collection:", [unit.var_set for unit in activeAtomicUnits])
    
    return activeAtomicUnits,inactiveAtomicUnits
    