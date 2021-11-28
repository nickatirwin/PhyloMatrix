# PhyloMatrix

Dependencies:
```
conda install -c conda-forge matplotlib
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c anaconda statsmodels
conda install -c anaconda scikit-learn

```
Usage:
```
from PhyloMatrix import PhyloMatrix
pm = PhyloMatrix()

# load clusters
pm.LoadMatrix(file='dump.out.reduced_eukaryotes_refproteomes.allVSall.dmnd.blastp.mci.I14',type='mcl')

# filter out clusters with less than n taxa
pm.Minimum_taxa(min_taxa=34)

# load cluster annotations
pm.Annotation.LoadAnnotations('COGannotations.tab')

### for association analysis
# run PCA and visualize
pm.Regression.PCA(binary=True)
pm.Regression.PCAPlot(PCx='PC1',PCy='PC2')

# assign the dependent variable (y)
pm.Regression.DependentVariable('phenotype.IFT74.tab')

# run the regression analysis
pm.Regression.Regression(a=0.05)
```
