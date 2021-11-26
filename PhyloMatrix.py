# to do
# 1. Compatibility for different orthogroup inputs?
# 2. Adjust PCA plot so that you can map different taxonomic groups

import sys
import pandas as pd
import numpy as np
from statsmodels.multivariate.pca import PCA
import statsmodels.api as sm
import matplotlib.pyplot as plt
from collections import Counter
from ete4 import NCBITaxa
from sklearn.preprocessing import StandardScaler

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
    
    def DependentVariable(self,dv=None,delim='\t',header=False):
        if not dv:
            raise ValueError('Provide a dependent variable file')
        try:
            if header == False:
                dv = open(dv,'r').readlines()
            else:
                dv = open(dv,'r').readlines()[1:]
            dv = {l.split(delim)[0]:int(l.split(delim)[1].strip()) for l in dv}
        except:
            raise ValueError("Cannot open dependent variable file: %s" % dv)
        PhyloMatrix.y = pd.DataFrame([dv[t] for t in PhyloMatrix.matrix.index],index=PhyloMatrix.matrix.index,columns=['y']) 
    
    #def Regression
    def Regression(self, a = 0.05, PC=None, remove_influentials=True):
        # select number of PCs to use: to 80% of variance explained
        if not PC:
            PC = 1
            while sum(PhyloMatrix.pca.eigenvals[0:PC])/sum(PhyloMatrix.pca.eigenvals) < 0.8:
                PC += 1
        PhyloMatrix.scale_matrix = pd.DataFrame(StandardScaler().fit_transform(PhyloMatrix.matrix), index=PhyloMatrix.matrix.index, columns=PhyloMatrix.matrix.columns)
        # run multiple regressions
        PhyloMatrix.regression_coefficients = {}
        PhyloMatrix.regression_pvalues = {}
        m = 0
        for cl in PhyloMatrix.matrix.columns:
            m += 1
            # progress bar
            sys.stdout.write('\r')
            sys.stdout.write(str(m)+'/'+str(len(PhyloMatrix.matrix.columns))+'     ')
            sys.stdout.flush()
            ####
            X = sm.add_constant(pd.concat([PhyloMatrix.scale_matrix[cl],PhyloMatrix.pca.values.iloc[:,:PC]],axis=1))
            # calculate initial model
            result = sm.GLM(PhyloMatrix.y, X).fit(cov='HC1')
            if remove_influentials == True:
                # identify and remove influential points
                influence = result.get_influence()
                n, k = X.shape[0], X.shape[1]
                cutoff_dffits = 2*(np.sqrt(k/n))
                X['dffits'] = influence.d_fittedvalues_scaled
                X['y'] = PhyloMatrix.y
                X_filt = X.loc[abs(X['dffits']) < cutoff_dffits]
                y_filt = X_filt['y']
                X_filt = X_filt.drop(columns=['y','dffits'])
                # recalculate model
                try:
                    result = sm.GLM(y_filt, X_filt).fit(cov='HC1')    
                except: # if the filtered results have perfect seperation, take previous analysis
                    result = sm.GLM(PhyloMatrix.y, X.drop(['dffits','y'],axis=1)).fit(cov='HC1')
            # record results
            PhyloMatrix.regression_coefficients[cl] = result.params[cl]
            PhyloMatrix.regression_pvalues[cl] = result.pvalues[cl]
        # return results
        a = a/len(PhyloMatrix.matrix.columns) # bonferonni correction
        regression_hits = [c for c in PhyloMatrix.matrix.columns if PhyloMatrix.regression_pvalues[c] < a]
        PhyloMatrix.regression_hits = pd.DataFrame(PhyloMatrix.matrix[regression_hits])
        PhyloMatrix.regression_hits.loc['coeff'] = [PhyloMatrix.regression_coefficients[c] for c in regression_hits]
        PhyloMatrix.regression_hits.loc['pvalue'] = [PhyloMatrix.regression_pvalues[c] for c in regression_hits]
        PhyloMatrix.regression_hits = PhyloMatrix.regression_hits.sort_values(by ='pvalue', axis=1)
    
    # plot results
    def RegressionPlot(self, a=0.05):
        a = a/len(PhyloMatrix.matrix.columns)
        unsig_p = [PhyloMatrix.regression_pvalues[p] for p in PhyloMatrix.regression_pvalues if PhyloMatrix.regression_pvalues[p] > a]
        unsig_c = [PhyloMatrix.regression_coefficients[p] for p in PhyloMatrix.regression_pvalues if PhyloMatrix.regression_pvalues[p] > a]
        sig_p = [PhyloMatrix.regression_pvalues[p] for p in PhyloMatrix.regression_pvalues if PhyloMatrix.regression_pvalues[p] < a]
        sig_c = [PhyloMatrix.regression_coefficients[p] for p in PhyloMatrix.regression_pvalues if PhyloMatrix.regression_pvalues[p] < a]
        plt.scatter(y=-1*(np.log10(unsig_p)),x=unsig_c,alpha=0.5)
        plt.scatter(y=-1*(np.log10(sig_p)),x=sig_c,alpha=0.5)
        plt.ylabel('-Log10(p value)')
        plt.xlabel('Coefficient')
        plt.show()