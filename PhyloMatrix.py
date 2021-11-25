# to do
# 1. Compatibility for different orthogroup inputs?
# 2. Adjust PCA plot so that you can map different taxonomic groups


import pandas as pd
import numpy as np
from statsmodels.multivariate.pca import PCA
import statsmodels.api as sm
import matplotlib.pyplot as plt
from collections import Counter
from ete4 import NCBITaxa

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
        "Generate a data matrix from a cluster file"
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
            if np.count_nonzero(PhyloMatrix.matrix[c]) < min_taxa: exclude.append(c)
        PhyloMatrix.matrix = PhyloMatrix.matrix.drop(exclude, axis=1)
    
    def PCA(self, binary=False):
        "Calculate a PCA using the gene matrix"
        # make the dataframe binary
        if binary == True:
            pca_matrix = pd.DataFrame(np.where(PhyloMatrix.matrix > 1, 1, 0),columns=PhyloMatrix.matrix.columns, index=PhyloMatrix.matrix.index).transpose()
        else:
            pca_matrix = PhyloMatrix.matrix.transpose()       
        # run the pca
        PhyloMatrix.pca = PCA(pca_matrix,standardize=True, method='eig')
        PhyloMatrix.pca.values = PhyloMatrix.pca.coeff.transpose()
        
    def PCAplot(self,PCx='PC1',PCy='PC2',annotate=False):
        "Plot for visualizing PCA"
        ncbi = NCBITaxa()
        # check taxonomic groups
        taxonomy, species = [], []
        for i in PhyloMatrix.pca.coeff.transpose().index:
            taxonomy.append(ncbi.translate_to_names([ncbi.get_lineage(int(i))[3]])[0])
            species.append(ncbi.translate_to_names([int(i)])[0])
        taxonomy_df = pd.DataFrame({"taxonomy":taxonomy,"species":species},index=PhyloMatrix.pca.coeff.columns)
        pca_fig = pd.concat([PhyloMatrix.pca.coeff.transpose(), taxonomy_df], axis = 1)
        # open figure and define axes
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1) 
        ax.set_xlabel(PCx+' ('+ str(round((PhyloMatrix.pca.eigenvals[int(PCx[-1])-1]/np.sum(PhyloMatrix.pca.eigenvals))*100,1)) + '% explained var.)')
        ax.set_ylabel(PCy+' ('+ str(round((PhyloMatrix.pca.eigenvals[int(PCy[-1])-1]/np.sum(PhyloMatrix.pca.eigenvals))*100,1)) + '% explained var.)')
        # make a list of taxonomic groups and colours
        groups = list(set(taxonomy))
        colours = (['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']*len(groups))[0:len(groups)]
        # add each taxonomic group to the plot
        for target, colour in zip(groups,colours):
            indicesToKeep = pca_fig['taxonomy'] == target
            ax.scatter(pca_fig.loc[indicesToKeep, pca_fig.columns[int(PCx[-1])]], pca_fig.loc[indicesToKeep, pca_fig.columns[int(PCy[-1])]], c=colour, s=50)
        # annotate points in axis
        if annotate == True:
            for idx, row in pca_fig.iterrows():
                ax.annotate(row['species'],(row[pca_fig.columns[int(PCx[-1])]], row[pca_fig.columns[int(PCy[-1])]]), fontsize=6 )
        ax.legend(groups)
        ax.grid()    
        # plot
        plt.show()
    
    #def AssignDependentVariable    
            
    #def PCFilter
    
    #def Regression