

class AtomicUnit:
    def __init__(self, var_set, surrogate_set_A, surrogate_set_B):
        self.var_set = var_set             #list, save the list of latent variable names in this atomic unit
        self.surrogate_set_A = surrogate_set_A     #len = len(self.var_set), list, save one measured surrogate set
        self.surrogate_set_B = surrogate_set_B     #len = len(self.var_set), list, save another measured surrogate set
        self.children = []            #list of list, save all children,include pure children and impure children
        self.overlapped_units = dict()    #dictionary, {atomic unit: number of overlapping}, save overlapped atomic units
        self.parents = []             #list of list, save all parents
    
    
    def add_child(self, child_set):
        if child_set not in self.children:
            self.children.append(child_set)    
   
    def update_children(self, child_set_list):
        for child_set in child_set_list:
            if child_set not in self.children:
                self.children.append(child_set)
    
    def remove_child(self, child_set):
        if child_set in self.children:
            self.children.remove(child_set) 
        
    def add_overlapped_units(self, atomicUnit, numOfOverlap):
        self.overlapped_units[atomicUnit] = numOfOverlap
    
    def remove_overlapped_units(self, atomicUnit):
        self.overlapped_units.pop(atomicUnit, 0)   
    
    def add_parent(self, parent_set):
        if parent_set not in self.parents:
            self.parents.append(parent_set)  
    
    def remove_parent(self, parent_set):
        if parent_set in self.parents:
            self.parents.remove(parent_set) 
    