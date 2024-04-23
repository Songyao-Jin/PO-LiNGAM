# [ICLR 2024] Structural Estimation of Partially Observed Linear Non-Gaussian Models: A Practical Approach with Identifiability
***
This is the official implementation of the paper [Structural Estimation of Partially Observed Linear Non-Gaussian Models: A Practical Approach with Identifiability](https://openreview.net/forum?id=nHkMm0ywWm), ICLR 2024.

### Overview
This project estimates the causal structure of the partial observed linear non-Gaussian acyclic model. The latent variable may be anywhere in the causal graph.


### Main_PO-LiNGAM Function
    Main_PO-LiNGAM.py : MainPOLiNGAM(data, alpha, maxNrOfParentUnits, maxAtomicUnitSize)
    Input:
        data: DataFrame type,  sample-by-dims size
        alpha: significance level of the independence test for three phases.
        maxNrOfParentUnits: the maximum possible number of parents of one atomic unit (set to to speed the inference of phase 1)
        maxAtomicUnitSize: the maximum possible number of variables in one atomic unit (set to to speed the inference of phase 2)
    Output:
        adj_matrix.txt : the adjacency matrix of the causal structure
        graph in plots file:  plot a causal graph 


### Test:
    One may modify the main() of Main_PO-LiNGAM.py to test our method.
    The default case is case 2. For simplicity, just run ``python Main_PO-LiNGAM.py``


### Notes
    Our method relies heavily on independence tests. 
    Here, HSIC-based Independence Test is used.

    Reference: Q. Zhang, S. Filippi, A. Gretton, and D. Sejdinovic, Large-Scale Kernel Methods for Independence Testing, Statistics and Computing, 2018.
