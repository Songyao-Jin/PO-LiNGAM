import pydot
from Atomic_unit import AtomicUnit
import time


def get_current_time():
    current_time = time.strftime('%Y-%m-%d_%H-%M-%S', time.localtime())
    return current_time

#draw Pyramid-shaped graph
def Make_mergedPyramid_graph(allAtomicUnits: list):
    graph = pydot.Dot(graph_type='digraph')

    for unit in allAtomicUnits:
        la = str()
        for varStr in unit.var_set:
            la = la + varStr +" "
        la = la[:-1]
        graph.add_node(pydot.Node(str(unit.var_set),label = la))
        

    for unit in allAtomicUnits:
        for parent in unit.parents:
            graph.add_edge(pydot.Edge(str(parent), str(unit.var_set)))
    

    added =[]   
    for unit in allAtomicUnits:
        for key, value in unit.overlapped_units.items():
            newset = {frozenset(key.var_set),frozenset(unit.var_set)}
            if newset not in added and key in allAtomicUnits:  
                graph.add_edge(pydot.Edge(str(unit.var_set), str(key.var_set), dir='none', label =value))
                added.append(newset) 
    

    current_time = get_current_time()
    graph.write_png(f"plots/merged_pyramid_graph_{current_time}.png")          
    graph.write_dot('plots/merged_pyramid_graph.dot')


#draw spring-shaped graph
def Make_springAndPyramid_graph(allAtomicUnits:list):

    graph = pydot.Dot(graph_type='digraph')
 
    clusterList = []

    for unit in allAtomicUnits:

        if len(unit.overlapped_units) != 0:
            cluster = pydot.Cluster(str(unit.var_set))
            graph.add_subgraph(cluster)
            clusterList.append(cluster)
            for varStr in unit.var_set:
                graph.add_node(pydot.Node(varStr))
                cluster.add_node(pydot.Node(varStr))
        else:
            for varStr in unit.var_set:
                graph.add_node(pydot.Node(varStr))

    for unit in allAtomicUnits:
        for parent in unit.parents:
            for var in unit.var_set:
                for paVar in parent:
                    if graph.get_node(paVar) !=[]: 
                        graph.add_edge(pydot.Edge(paVar, var))
    
    added =[]   
    for unit in allAtomicUnits:
        for key, value in unit.overlapped_units.items():
            newset = {frozenset(key.var_set),frozenset(unit.var_set)}
            if newset not in added and key in allAtomicUnits:   
                clu_unit=clu_key=pydot.Cluster("dummy")
                for clu in clusterList:
                    name = clu.get_name()  
                    if name == '"'+'cluster_'+str(unit.var_set)+'"':
                        clu_unit = clu
                    if name == '"'+'cluster_'+str(key.var_set)+'"':
                        clu_key = clu   
                if clu_key.get_name() != '"cluster_dummy"':    
                    graph.add_edge(pydot.Edge(clu_unit, clu_key, color="blue", dir='none', label =value))
                    added.append(newset) 
    

    graph.set_layout('dot')
    graph.write_png("plots/pyramid_graph.png")      
    graph.write_dot('plots/pyramid_graph.dot') 
    
    graph.set_layout('fdp')
    graph.write_png("plots/spring_graph.png")      
    graph.write_dot('plots/spring_graph.dot')     




