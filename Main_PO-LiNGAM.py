import numpy as np
import pandas as pd
import warnings
import MakeGraph 

from Atomic_unit import AtomicUnit
import Phase1_IdentifyCausalStructureAmongKnownAtomicUnits as Phase1
import Phase2_FindNewAtomicVariableSetByClustering as Phase2
import Phase3_CheckOverlapBetweenAtomicUnits as Phase3
import SimulationData as SD #simulate data


def MainPOLiNGAM(data, alpha=[0.05,0.01,0.01], maxNrOfParentUnits = 4, maxAtomicUnitSize = 4):     #set a contraints of maxAtomicUnitSize to reduce calculation complexity
    
    indexs=list(data.columns)
    activeAtomicUnits = []
    inactiveAtomicUnits = []
    activeAtomicUnits_copy = activeAtomicUnits.copy()
    
    
    for col in data.columns:
        instance = AtomicUnit([col],[col],[col])    
        activeAtomicUnits.append(instance)
    
    
    currentNumOfLatentVariables = 0
    numOfIter = 0
    while len(activeAtomicUnits) > 1 and activeAtomicUnits_copy != activeAtomicUnits:    #若是在聚类的时候还是等于2，那没办法识别是否他们会归于一个 #also # #Check if code need to exit the while loop
        
        activeAtomicUnits_copy = activeAtomicUnits.copy()
        numOfIter +=1
        # phase1
        print("++++++++ Phase 1 start +++++++")
        activeAtomicUnits, inactiveAtomicUnits = Phase1.IdentifyCausalStructureAmongKnownAtomicUnits(data, activeAtomicUnits, inactiveAtomicUnits, alpha[0], maxNrOfParentUnits)
        
        # phase2
        print('')
        print("++++++++ Phase 2 start +++++++")
        newGeneratedAtomicUnits, activeAtomicUnits, inactiveAtomicUnits, currentNumOfLatentVariables  = Phase2.FindNewAtomicVariableSetByClustering(data, activeAtomicUnits, inactiveAtomicUnits, currentNumOfLatentVariables, alpha[1], maxAtomicUnitSize)

        # phase3
        print('')
        print("++++++++ Phase 3 start +++++++")
        activeAtomicUnits  = Phase3.CheckOverlapBetweenAtomicUnits(data, activeAtomicUnits, inactiveAtomicUnits, newGeneratedAtomicUnits, alpha[2])
        print('')
        print("++++++++++++ Summary ++++++++++++++++")
        print("Current Atomic Unit Collection:", end='')
        for act in activeAtomicUnits:
            print(act.var_set, end='')
        print("")
        print("++++++++++++ The",str(numOfIter)+"th","iteration end ++++++++++++++++")
        print("\n")
        
    
    
    direct_causal_relation ={}
    for unit in inactiveAtomicUnits:
        for parents in unit.parents:
            for pa in parents:
                if pa not in direct_causal_relation.keys():
                    direct_causal_relation[pa] = unit.var_set
                else:
                    direct_causal_relation[pa] = direct_causal_relation[pa] + unit.var_set   

    m=UpdateGraph(indexs,direct_causal_relation)
    file_name = 'adj_matrix.txt'
    with open(file_name, 'a') as f:
        df_str = m.to_string(index=True)
        f.write(df_str + '\n\n') 
    



    allAtomicUnits = inactiveAtomicUnits + activeAtomicUnits
    
    MakeGraph.Make_mergedPyramid_graph(allAtomicUnits)
    
    print("Causal relations:")
    for unit in allAtomicUnits:
        print("The atomic unit",unit.var_set,"has children:",unit.children)
    
    return None
    


def UpdateGraph(obserd,LatentIndex):
    key =LatentIndex.keys()
    Variables=[]
    for i in obserd:
        Variables.append(i)
    for i in key:
        if i not in Variables:
            Variables.append(i)
    n=len(Variables)
    indexs=Variables
    matrix=pd.DataFrame(np.zeros((n,n),dtype=np.int32))
    matrix.columns=indexs
    matrix.index=indexs
    for i in key:
        clu=LatentIndex[i]
        for j in clu:
            matrix[i][j]=1
    print(matrix)
    return matrix






#simulation: Randomly generate simulation data to validate our method
def main():

    for i in range(1):
        #generate simulation data
        data=SD.CaseS2(10000)
        #set the alpha value
        alpha_phase1 = 0.002
        alpha_phase2 = 0.002
        alpha_phase3 = 0.002
        alpha = [alpha_phase1,alpha_phase2,alpha_phase3]
        #set contraints to reduce calculation complexity
        maxAtomicUnitSize = 2
        maxNrOfParentUnits = 2
        #esimate causal structure from observed data, plot the graph
        MainPOLiNGAM(data, alpha, maxNrOfParentUnits, maxAtomicUnitSize)
        

if __name__ == '__main__':
    main()