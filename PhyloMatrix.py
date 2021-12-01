# to do
# 1. Compatibility for different orthogroup inputs?
# 2. Adjust PCA plot so that you can map different taxonomic groups
# 3. Speed up annotation subsampling using regex: pm.annotation_matrix.go.str.contains('^go1')

import sys
import pandas as pd
import numpy as np
from statsmodels.multivariate.pca import PCA
import statsmodels.api as sm
import matplotlib.pyplot as plt

import seaborn

import mpld3
from mpld3 import plugins

from collections import Counter
from sklearn.preprocessing import StandardScaler

from ete4 import NodeStyle
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace
from ete4.treeview import random_color


class PhyloMatrix(object):
    def __init__(self, cluster_file=None):
        # initialize objects
        PhyloMatrix.cluster_file = None
        PhyloMatrix.matrix = None
        PhyloMatrix.tree = None
        PhyloMatrix.annotation_matrix = pd.DataFrame()
    
    # MATRIX PARSER
    def LoadMatrix(self, file=None, type='table', min_taxa=1, header=0, sep='\t'):
        '''Load a matrix from clusters'''

        def Minimum_taxa(d,families=None,min_taxa=min_taxa):    
            # filter clusters with not enough taxa
            exclude = []
            for m in range(0,len(d[next(iter(d))])):
                taxa = 0
                for t in d:
                    if d[t][m] > 0:
                        taxa += 1
                if taxa < min_taxa:
                    exclude.append(m)
            for t in d:
                d[t] = list(np.delete(np.array(d[t]),exclude))
            if families:
                families = list(np.delete(np.array(families),exclude))
            return d,families

        # check for input file
        if file:
            try:
                PhyloMatrix.cluster_file = open(file,'r').read()
            except:
                raise ValueError("Cannot open cluster file: %s" % file)
        else:
            raise ValueError("Please provide cluster file")

        # check input file type and load    
        if type.lower() == 'mcl': # load in MCL output
            # record all taxa
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
            taxa_d = Minimum_taxa(taxa_d,min_taxa=min_taxa)
            PhyloMatrix.matrix = pd.DataFrame.from_dict(taxa_d, columns=['C'+(6-len(str(x)))*'0'+str(x) for x in range(1,n+1)],orient='index')   
        
        elif type.lower() == 'eggnog': # load in eggnog family members
            # record all taxa
            taxa = []
            for l in PhyloMatrix.cluster_file.split('\n')[:-1]:
                taxa += [t.split('.')[0] for t in l.split('\t')[4].split(',')]
            taxa = set(taxa)
            # make a dictionary
            n = len(PhyloMatrix.cluster_file.split('\n')[:-1])
            taxa_d = {t:[0]*n for t in taxa}
            # count gene frequencies
            m, family = 0, []
            for l in PhyloMatrix.cluster_file.split('\n')[:-1]:
                counts = Counter([t.split('.')[0] for t in l.split('\t')[4].split(',')])
                for t in counts:
                    taxa_d[t][m] = counts[t]
                family.append(l.split('\t')[1])
                m += 1
            taxa_d,family = Minimum_taxa(taxa_d,families=family,min_taxa=min_taxa)
            PhyloMatrix.matrix = pd.DataFrame.from_dict(taxa_d, columns=family,orient='index')  

        elif type.lower() == 'table': # load in table as a matrix
            # read in the dataframe using pandas
            PhyloMatrix.matrix = pd.read_table(file,header=header, index_col=0, sep=sep)
            # if there are no headers, provide them
            if header != 0:
                PhyloMatrix.matrix.columns = ['C'+(6-len(str(x)))*'0'+str(x) for x in range(1,len(PhyloMatrix.cluster_file.split('\n')[:-1])+1)]

    class PCA(object):
        def __init__(self, binary=True, standardize=True, method='eig', normalize=True):
            "Calculate a PCA using the gene matrix"
            # make the dataframe binary
            if binary == True:
                pca_matrix = pd.DataFrame(np.where(PhyloMatrix.matrix > 1, 1, 0),columns=PhyloMatrix.matrix.columns, index=PhyloMatrix.matrix.index).transpose()
            else:
                pca_matrix = PhyloMatrix.matrix.transpose()       
            # run the pca - max 100 components
            if len(pca_matrix.columns) > 100:
                n = 100
            else:
                n = len(pca_matrix.columns)
            PhyloMatrix.pca = PCA(pca_matrix,standardize=standardize, method=method, normalize=normalize, ncomp=n)
            PhyloMatrix.pca.values = PhyloMatrix.pca.coeff.transpose()

        def plot(PCx='PC1',PCy='PC2',annotate=False):
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
            #mpld3.save_html(fig, "pca.html")
            plt.show()


    class ProfileTree(object):
        def __init__(self):
            pass
        
        def plotProfile(self):
            def get_profile_layout():
                def get_color(cmap, value):
                    cmap = cm.get_cmap(cmap)
                    return cmap(value)
                def layout_fn(node):
                    # iterate through profile
                    profile = node.props.get('profile')
                    #if profile:
                    for i, value in enumerate(profile):
                        col = rgb2hex(get_color("binary", value))

                        w = 50
                        h = 50

                        f = RectFace(w, h, color=col, padding_x=1, padding_y=1)

                        if node.is_leaf():
                            node.add_face(f,  column=2+i, position="aligned")
                        else:

                            node.add_face(f,  column=2+i, position="aligned",  collapsed_only=True)

                layout_fn.__name__ == "Profile Matrix"
    
                layout_fn.contains_aligned_face = True
                
                return layout_fn

            def get_pheno_layout():
                def get_color(cmap, value):
                    cmap = cm.get_cmap(cmap)
                    return cmap(value)

                def layout_fn(node):
                    if node.props.get('pheno'):
                        try:
                            i = 2 + len(node.props.get('profile'))
                        except:
                            i = 2
                        col = rgb2hex(get_color("Reds", node.props.get('pheno') * 256))
                        f = RectFace(50, 50, color=col, padding_x=10, padding_y=1)

                    if node.is_leaf():
                        node.add_face(f,  column=i, position="aligned")
                    else:
                        node.add_face(f,  column=i, position="aligned",  collapsed_only=True)

                layout_fn.__name__ = 'Phenotype'
                layout_fn.contains_aligned_face = True
                return layout_fn


            ncbi = NCBITaxa()
            species = matrix.index

            t = ncbi.get_topology(taxids)
            t.convert_to_ultrametric()

            # Annotate tree
            for node in t.traverse():
                if node.is_leaf():
                    node.add_prop("profile", list(PhyloMatrix.matrix.loc[node.name]))
                    # treating phenotype as number for now
                    node.add_prop("pheno", int(PhyloMatrix.y.loc[node.name]))
                    
            else:
                node.add_prop("profile", list(PhyloMatrix.matrix.loc[node.get_leaf_names()].mean(axis=0)))
                node.add_prop("pheno", float(PhyloMatrix.y.loc[node.get_leaf_names()].mean()))
            
            ts = TreeStyle()
            ts.show_branch_length = False
            ts.show_branch_support = False
            ts.show_leaf_name = False

            # NEW: dd layouts accessible from front end
            profile_layout = get_profile_layout()
            profile_layout.__name__ = 'Profile Matrix'
            
            ts.layout_fn = [profile_layout, get_pheno_layout]
            
            t.explore(tree_name="test", tree_style=ts)



            
    class Regression(object):
        def __init__(self):
            pass

        def DependentVariable(dv=None,delim='\t',header=False):
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
        def Regression(a = 0.05, PC=None, remove_influentials=True):
            a = a/len(PhyloMatrix.matrix.columns) # bonferonni correction
            # select number of PCs to use: to 70% of variance explained
            if not PC:
                PC = 1
                while sum(PhyloMatrix.pca.eigenvals[0:PC])/sum(PhyloMatrix.pca.eigenvals) < 0.7:
                    PC += 1
            PhyloMatrix.scale_matrix = pd.DataFrame(StandardScaler().fit_transform(PhyloMatrix.matrix), index=PhyloMatrix.matrix.index, columns=PhyloMatrix.matrix.columns)
            # run multiple regressions
            PhyloMatrix.regression_coefficients = {}
            PhyloMatrix.regression_pvalues = {}
            PhyloMatrix.regression_association = {}
            m = 0
            for cl in PhyloMatrix.matrix.columns:
                m += 1
                # progress bar
                sys.stdout.write('\r')
                sys.stdout.write(str(m)+'/'+str(len(PhyloMatrix.matrix.columns))+'     ')
                sys.stdout.flush()
                ####
                try:
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
                    if (PhyloMatrix.regression_coefficients[cl] > 0) and (PhyloMatrix.regression_pvalues[cl] < a):
                        PhyloMatrix.regression_association[cl] = 'positive'
                    elif (PhyloMatrix.regression_coefficients[cl] < 0) and (PhyloMatrix.regression_pvalues[cl] < a):
                        PhyloMatrix.regression_association[cl] = 'negative'
                    else:
                        PhyloMatrix.regression_association[cl] = 'NA'
                except:
                    print('Error perfect correlation: %s'% cl)
            # return results
            correlated_clusters = {c:[PhyloMatrix.regression_coefficients[c],PhyloMatrix.regression_pvalues[c],PhyloMatrix.regression_association[c]] for c in PhyloMatrix.regression_coefficients.keys()}
            PhyloMatrix.correlated_clusters = pd.DataFrame.from_dict(correlated_clusters,orient='index',columns=['coeff','pvalue','association'])
            # add the regression output to the annotation matrix
            if PhyloMatrix.annotation_matrix.empty:
                PhyloMatrix.annotation_matrix = correlated_clusters
            else:
                PhyloMatrix.annotation_matrix = PhyloMatrix.annotation_matrix.merge(PhyloMatrix.correlated_clusters, left_index=True, right_index=True)

        # plot results
        def RegressionPlot(a=0.05):

            seaborn.relplot(data=PhyloMatrix.correlated_clusters,x='coeff',y=-np.log(PhyloMatrix.correlated_clusters.pvalue),hue='association')
            '''    
            css = """
            table
            {
            border-collapse: collapse;
            }
            th
            {
            color: #ffffff;
            background-color: #000000;
            }
            td
            {
            background-color: #cccccc;
            }
            table, th, td
            {
            font-family:Arial, Helvetica, sans-serif;
            border: 1px solid black;
            text-align: right;
            }
            """
            
            a = a/len(PhyloMatrix.matrix.columns)
            unsig_p = [PhyloMatrix.regression_pvalues[p] \
                       for p in PhyloMatrix.regression_pvalues \
                       if PhyloMatrix.regression_pvalues[p] > a]
            unsig_c = [PhyloMatrix.regression_coefficients[p] \
                       for p in PhyloMatrix.regression_pvalues \
                       if PhyloMatrix.regression_pvalues[p] > a]
            sig_p = [PhyloMatrix.regression_pvalues[p] \
                     for p in PhyloMatrix.regression_pvalues \
                     if PhyloMatrix.regression_pvalues[p] < a]
            sig_c = [PhyloMatrix.regression_coefficients[p] \
                     for p in PhyloMatrix.regression_pvalues \
                     if PhyloMatrix.regression_pvalues[p] < a]


            fig = plt.figure(figsize = (8,8))
            
            x = unsig_c
            y = -1*(np.log10(unsig_p))            
            points1 = plt.scatter(y=x,x=y,alpha=0.5)
            labels1 = len(x) * ['test']

            x = sig_c
            y = -1*(np.log10(sig_p))
            points2 = plt.scatter(y=y,x=x,alpha=0.5)
            labels2 = len(x) * ['test']
            plt.ylabel('-Log10(p value)')
            plt.xlabel('Coefficient')            

            """
            tooltip1 = plugins.PointHTMLTooltip(points1[0],
                                                labels,
                                                voffset=10,
                                                hoffset=10,
                                                css=css)

            tooltip2 = plugins.PointHTMLTooltip(points2[0],
                                                labels,
                                                voffset=10,
                                                hoffset=10,
                                                css=css)
            
            plugins.connect(fig, tooltip1)
            plugins.connect(fig, tooltip2)
            """
            mpld3.save_html(fig, "regression.html")
            '''
            plt.show()

    def LoadAnnotations(self, file=None, sep='\t', header=True):
        annot = open(file,'r').readlines()
        if header == True:
            a_d = {a.split(sep)[0]:a.strip().split('\t')[1:] for a in annot[1:]}
        else:
            a_d = {a.split(sep)[0]:a.strip().split('\t')[1:] for a in annot}
        PhyloMatrix.annotation_matrix = pd.DataFrame.from_dict(a_d,orient='index',columns=annot[0].strip().split('\t'))
    
    def Subset(self, targets=None,features=None):
        '''provide a list of targets and a list of features - search for the targets in each feature'''
        if type(targets) == str: targets = [targets]
        if type(features) == str: features = [features]
        sub = dict(zip(features,targets))
        clusters = []
        for s in sub:
            clusters.append(PhyloMatrix.annotation_matrix.loc[PhyloMatrix.annotation_matrix[s].isin([sub[s]])].index.to_list())
        clusters = set.intersection(*map(set,clusters))
        PhyloMatrix.extracted_annotations = PhyloMatrix.annotation_matrix.loc[clusters,:] # Just for testing
        return PhyloMatrix.matrix.loc[:,clusters]

