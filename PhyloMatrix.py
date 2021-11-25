# to do
# 1. Compatibility for different orthogroup inputs?



import pandas as pd
import numpy as np
from statsmodels.multivariate.pca import PCA
import statsmodels.api as sm
import matplotlib.pyplot as plt
from collections import Counter

class PhyloMatrix(object):
    def __init__(self, cluster_file=None):
        if cluster_file:
            try:
                PhyloMatrix.cluster_file = open(cluster_file,'r').read()
            except:
                raise ValueError("Cannot open cluster file: %s" % cluster_file)
        if not cluster_file:
            raise ValueError("Please provide cluster file")
    
    def GenerateMatrix(self, cluster_file=None, min_taxa=None):
        # record all samples
        taxa = set([t.split('.')[0] for t in PhyloMatrix.cluster_file.split('\t')])
        n = len(PhyloMatrix.cluster_file.split('\n')[:-1])
        taxa_d = {t:[0]*n for t in taxa}
        # count gene frequencies
        m = 0
        for c in PhyloMatrix.cluster_file.split('\n')[:-1]:
            counts = Counter([s.split('.')[0] for s in c.split('\t')])
            for t in counts:
                taxa_d[t][m] = counts[t] 
            m += 1
        # make dataframe
        PhyloMatrix.matrix = pd.DataFrame(taxa_d, index=['C'+(6-len(str(x)))*'0'+str(x) for x in range(1,n+1)]).transpose()
        # filter cols with not enough taxa
        if not min_taxa:
            min_taxa = 0
        exclude = []
        for c in PhyloMatrix.matrix.columns:
            if np.count_nonzero(PhyloMatrix.matrix[c]) < minimum_taxa: exclude.append(c)
        PhyloMatrix.matrix = PhyloMatrix.matrix.drop(exclude, axis=1)
    
    def PCA(self):
        if not PhyloMatrix.matrix:
            raise ValueError("No matrix found: Run GenerateMatrix()")
        PhyloMatrix.pca = PCA(PhyloMatrix.matrix,standardize=True)
    
    #def PCFilter
    
    #def AssignDependentVariable
    
    #def Regression