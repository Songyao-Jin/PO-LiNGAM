import numpy as np
import pandas as pd
from itertools import combinations
import itertools
import networkx as nx





def bfs(graph, start):
    visited, queue = set(), [start]
    while queue:
        vertex = queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            queue.extend(graph[vertex] - visited)
    return visited

def connected_components(G):
    seen = set()
    for v in G:
        if v not in seen:
            c = set(bfs(G, v))
            yield c
            seen.update(c)

def graph(edge_list):
    result = {}
    for source, target in edge_list:
        result.setdefault(source, set()).add(target)
        result.setdefault(target, set()).add(source)
    return result



def merge_list_modified(list_of_lists):
    merged_lists = []
    
    for sublist in list_of_lists:
        merged = False
        for existing_list in merged_lists:
            if any(item in existing_list for item in sublist):
                existing_list.extend(item for item in sublist if item not in existing_list)
                merged = True
                break
        if not merged:
            merged_lists.append(sublist)
    
    return merged_lists

# # Example usage:
# list_of_lists = [['A', 'B', 'C'], ['C', 'D', 'E'], ['F', 'G']]
# merged_lists = merge_list_modified(list_of_lists)
# print(merged_lists)






def merge_list(L2):
    l=L2.copy()

    edges = []
    s = list(map(set, l))
    for i, j in combinations(range(len(s)), r=2):

        if s[i].intersection(s[j]):
            edges.append((i, j))

    G = graph(edges)

    result = []
    unassigned = list(range(len(s)))

    for component in connected_components(G):

        union = set().union(*(s[i] for i in component))

        result.append(sorted(union))

        unassigned = [i for i in unassigned if i not in component]


    result.extend(map(sorted, (s[i] for i in unassigned)))

    return result


#merge withour change any order
def merge_listNew(L2):
    
    l=L2.copy()

    edges = []
    s = list(map(set, l))
    for i, j in combinations(range(len(s)), r=2):
        if s[i].intersection(s[j]):
            edges.append((i, j))

    G = graph(edges)

    result = []
    unassigned = list(range(len(s)))

    for component in connected_components(G):

        component=list(component)

        Cluster=l[component[0]]

        adj= s[0].intersection(s[1])
        for i in component:
            if i == 0:
                continue
            temp=l[i]
            for j in list(adj):
                if j in temp:
                    temp.remove(j)

            for j in range(0,len(temp)):
                if j < (len(temp)/2):
                    Cluster.insert(0,temp[j])
                else:
                    Cluster.append(temp[j])

        result.append(Cluster)



        unassigned = [i for i in unassigned if i not in component]


    result.extend(map(sorted, (s[i] for i in unassigned)))

    return result




def main():
    L=[['x1','x2'],['x3','x4'],['x5','x6','x2'],['x1','x9']]
    print(merge_list(L))


if __name__ == '__main__':
    main()
