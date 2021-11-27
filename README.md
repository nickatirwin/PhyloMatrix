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
pm.LoadMatrix(file='dump.out.reduced_eukaryotes_refproteomes.allVSall.dmnd.blastp.mci.I14',type='mcl')
pm.Minimum_taxa(34)
pm.Annotation.LoadAnnotations('COGannotations.tab')
pm.Regression.PCA(binary=True)
pm.Regression.DependentVariable('phenotype.IFT74.tab')
pm.Regression.Regression()
