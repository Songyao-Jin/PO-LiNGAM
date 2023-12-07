import itertools
import GIN2
import Utils_V2 as Utils 
from Atomic_unit import AtomicUnit


#return the combination of list
def FindCombination(Lists,N):
    return itertools.combinations(Lists,N)


#GIN test by fast HSIC
#X and Z are list, e.g., X=['X1','X2'] Z=['X3']
#Data.type=Pandas.DataFrame, where data.columns=['x1','x2',...]
def GIN(X,Z,data,alpha=0.05):
    return GIN2.GIN(X,Z,data,alpha)



def FindNewAtomicVariableSetByClustering(data, activeAtomicUnits, inactiveAtomicUnits, currentNumOfLatentVariables, alpha, maxAtomicUnitSize):

    newGeneratedAtomicUnits = []
    haveChildUnitCollection= set()    
    
    if len(activeAtomicUnits) <=1:
        return newGeneratedAtomicUnits, activeAtomicUnits, inactiveAtomicUnits, currentNumOfLatentVariables 
            
    
    # one-factor cluster
    Grlen=2         
    pairSets=[]    
    pureSets=[]
    tempUnknownPairSets=[]
    #1-factor : pair-with-pair learning for pair set (pure or impure set)
    Sets_P=FindCombination(activeAtomicUnits,Grlen)
    for P in Sets_P:
        Z_set = activeAtomicUnits.copy()
        Y_set_surrogate = list()
        Z_set_surrogate = list()
        
        
        if P[0] in P[1].overlapped_units or P[1] in P[0].overlapped_units:
            continue
        
        for y in P:
            Z_set.remove(y)   
            for key, value in y.overlapped_units.items():
                if len(key.var_set) == value and value <= len(y.var_set):   
                    if key in Z_set:
                        Z_set.remove(key)
            Y_set_surrogate.append(y.surrogate_set_A[0])
        
        for z in Z_set:
            Z_set_surrogate.append(z.surrogate_set_A[0])

        if GIN(Y_set_surrogate,Z_set_surrogate,data,alpha):    
            pairSets.append(list(P))
            print([unit.var_set for unit in P], [unit.var_set for unit in Z_set], "satisfies GIN")



    temp=Utils.merge_list_modified(pairSets)
    pairSets=temp 
    for se in pairSets:
        if len(se)>=3:
            pureSets.append(se)
        else:
            tempUnknownPairSets.append(se)    
            
            

    purePairSets, impurePairSets = IdentifyUnknownPairSet(data, tempUnknownPairSets, activeAtomicUnits, alpha=0.01)
    pureSets= pureSets + purePairSets


    for Pset in pureSets:
        fullPureChildren = Pset    
        pureChildsetA = set([Pset[0]])
        surrogate_pureChildSetB = [Pset[1].surrogate_set_A[0]]
        currentNumOfLatentVariables, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits = CreateNewAtomicVariableSet(Grlen-1, currentNumOfLatentVariables, fullPureChildren, pureChildsetA, surrogate_pureChildSetB, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits)
     
    
    
    # multi-factor cluster
    Grlen += 1     
    clusterList = []   
    while(Grlen-1 <= maxAtomicUnitSize ):  
        surrogate_activeAtomicUnits = []
        for unit in activeAtomicUnits:
            surrogate_activeAtomicUnits.extend(unit.surrogate_set_A)
        
        # multi-factor, 3 or more
        Sets_P=FindCombination(surrogate_activeAtomicUnits,Grlen) 
        for Y_set_surrogate in Sets_P:
            Y_set_surrogate = list(Y_set_surrogate)
            units=set()     
            haveChildUnit = False
            for act in activeAtomicUnits:
                for var in Y_set_surrogate:
                    if var in act.surrogate_set_A: 
                        units.add(act) 
                        if act in haveChildUnitCollection:       
                            haveChildUnit = True
            if haveChildUnit:
                continue   
               
            
            if len(units)>1 and units not in clusterList:
                Z_set = activeAtomicUnits.copy()
                Z_set_surrogate = list()
                for at in list(units):     
                    try:
                        Z_set.remove(at)
                    except ValueError:
                        print("The unit",at.var_set,"may missed from Z_set.")
                        print("The original Z_set is,",[unit.var_set for unit in activeAtomicUnits])
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
                    haveChildUnit,  haveChildUnitCollection = CheckWhetherHavingPartialParent(units, Y_set_surrogate, Z_set_surrogate, haveChildUnitCollection,data, alpha)
                    if not haveChildUnit:
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
                                # print("difference_unit: ",difference_unit)
                                if len(cluster)==len(baseCluster) and len(difference_unit) == 1:   
                                    pureChildsetA.update(difference_unit)
                                    # print("cluster: ", [x.var_set for x in cluster])
                        
                        if len(pureChildsetA) >= (Grlen-1): 
                            # print("pureChildren: ", [x.var_set for x in pureChild])
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
                                print("mergedClusterListTillNow: ",len(mergedClusterList))
                                
                                
                                settled = []; unsettled = []
                                for child in fullPureChildren:
                                    if child not in activeAtomicUnits:
                                        settled.append(child)
                                    else:
                                        unsettled.append(child)
                                if len(settled) == 0:
                                    currentNumOfLatentVariables, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits = CreateNewAtomicVariableSet(Grlen-1, currentNumOfLatentVariables, fullPureChildren, pureChildsetA, surrogate_pureChildSetB, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits)            
                                else:  
                                    for unit in  newGeneratedAtomicUnits:      
                                        if unit.var_set == settled[0].parents[0]:
                                            for un in unsettled:
                                                unit.add_child(un.var_set)
                                                un.add_parent(unit.var_set)
                                                activeAtomicUnits.remove(un)
                                                inactiveAtomicUnits.append(un)
                                            break   
                                break
        
        Grlen += 1        
    print('\n')   
    
    activeAtomicUnits.extend(newGeneratedAtomicUnits)
    
    print("Phase2 End:")
    print("Current Active Atomic Unit Collection: ")
    for unit in activeAtomicUnits:
        print(" unit :", unit.var_set)
    print("Current Inactive Atomic Unit Collection: ")
    for unit in inactiveAtomicUnits:
        print(" unit :", unit.var_set)
    print("")

    return newGeneratedAtomicUnits, activeAtomicUnits, inactiveAtomicUnits, currentNumOfLatentVariables                       
                


        
        







def CheckWhetherHavingPartialParent(unitSet, Y_set_surrogate, Z_set_surrogate, haveChildUnitCollection, data, alpha):
    
    haveChildUnit = False
    for unit in unitSet:
        new_Z_set_surrogate = Z_set_surrogate + unit.surrogate_set_B   
        if GIN(Y_set_surrogate,new_Z_set_surrogate,data,alpha):
            haveChildUnitCollection.add(unit)
            haveChildUnit = True
        
    return haveChildUnit,  haveChildUnitCollection


             
              


def CreateNewAtomicVariableSet(varSetSize, currentNumOfLatentVariables,fullPureChildren, pureChildsetA, surrogate_pureChildSetB, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits):
    
    var_set = ['L'+str(c) for c in range(currentNumOfLatentVariables+1,currentNumOfLatentVariables+varSetSize+1)]
    
    surrogate_set_A = []
    pureChildsetA_list = list(pureChildsetA)
    for unit in pureChildsetA_list[0:varSetSize]:
        surrogate_set_A.append(unit.surrogate_set_A[0])
    
    surrogate_set_B = surrogate_pureChildSetB[0:varSetSize]
        
    newAtomicUnit = AtomicUnit(var_set,surrogate_set_A,surrogate_set_B)    
    

    newGeneratedAtomicUnits.append(newAtomicUnit)
    for unit in fullPureChildren:
        unit.add_parent(newAtomicUnit.var_set) 
        newAtomicUnit.add_child(unit.var_set)

        if unit in activeAtomicUnits:
            activeAtomicUnits.remove(unit)
        else:
            print(f"{unit} is not in the list.")

        inactiveAtomicUnits.append(unit)
    

    print("NewAtomicVariableSet :", newAtomicUnit.var_set)
    print("Current Active Atomic Unit Collection: ")
    currentAll = activeAtomicUnits + [newAtomicUnit]
    for unit in currentAll:
        print(" unit :", unit.var_set)
    print("Current Inactive Atomic Unit Collection: ")
    for unit in inactiveAtomicUnits:
        print(" unit :", unit.var_set)
    print("")
    return currentNumOfLatentVariables+varSetSize, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits      
        
    
   
    
def IdentifyUnknownPairSet(data, tempUnknownPairSets, activeAtomicUnits, alpha=0.01):
    purePairSets=[]
    impurePairSets=[]
    for un_set in tempUnknownPairSets:
        un_set_surrogate_A = []
        un_set_surrogate_B = []
        restActiveAtomicUnits = activeAtomicUnits.copy()
        for u in un_set:
            un_set_surrogate_A.append(u.surrogate_set_A[0])
            un_set_surrogate_B.append(u.surrogate_set_B[0]) 
            restActiveAtomicUnits.remove(u) 
          
        restActiveAtomicUnits_surrogate = []
        for res in restActiveAtomicUnits:
            restActiveAtomicUnits_surrogate.append(res.surrogate_set_A[0])
            
                      
        impureFlag = False
        pairs =FindCombination(restActiveAtomicUnits_surrogate,2)
        for pair in pairs:
            if (GIN([un_set_surrogate_A[0], un_set_surrogate_A[1], pair[0]],[un_set_surrogate_B[0], pair[1]], data, alpha) ^ GIN([un_set_surrogate_A[0], un_set_surrogate_A[1], pair[0]],[un_set_surrogate_B[1], pair[1]], data, alpha)):  
                if (GIN([un_set_surrogate_B[0], un_set_surrogate_B[1], pair[0]],[un_set_surrogate_A[0], pair[1]], data, alpha) ^ GIN([un_set_surrogate_B[0], un_set_surrogate_B[1], pair[0]],[un_set_surrogate_A[1], pair[1]], data, alpha)) and (GIN([un_set_surrogate_B[0], un_set_surrogate_B[1], pair[1]],[un_set_surrogate_A[0], pair[0]], data, alpha) ^ GIN([un_set_surrogate_B[0], un_set_surrogate_B[1], pair[1]],[un_set_surrogate_A[1], pair[0]], data, alpha))  and  (GIN([un_set_surrogate_A[0], un_set_surrogate_A[1], pair[1]],[un_set_surrogate_B[0], pair[0]], data, alpha) ^ GIN([un_set_surrogate_A[0], un_set_surrogate_A[1], pair[1]],[un_set_surrogate_B[1], pair[0]], data, alpha)):              
                    print("impure set: "+ str([un_set[0].var_set,un_set[1].var_set])+" with "+str(pair))
                    impureFlag = True
                    break
        if impureFlag == True:
            impurePairSets.append(un_set)
        else:
            purePairSets.append(un_set)
            
    return purePairSets, impurePairSets
        
    
            
        
            
        
    