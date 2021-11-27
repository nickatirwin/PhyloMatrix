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

# load in the cluster file
phylo = PhyloMatrix('dump.out.reduced_eukaryotes_refproteomes.allVSall.dmnd.blastp.mci.I14')

# generate a cluster matrix
phylo.GenerateMatrix(min_taxa=34)

# calculate PCA and visualize
phylo.PCA(binary=True)
phylo.PCAplot(PCx='PC1',PCy='PC2',annotate=False)

# define the dependent variable (phenotype)
phylo.DependentVariable('phenotype.IFT74.tab')

# run a regression analysis to identify correlating genes and visualize
phylo.Regression()
phylo.RegressionPlot()
