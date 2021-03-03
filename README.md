# PhyloEpiGenomics

This package provides basic phylogenetic tree reconstruction algorithms for genomic or epigenomic data such as maximum likelihood, distance-based methods and parsimony. The library also allows to simulate sequence evolution based on evolutionary models and basic phylogenetic hypothesis testing.

## Installation

```R
install.packages("devtools", repos="http://cloud.r-project.org", clean=T)
devtools::install_github("hoffmann-lab/PhyloEpiGenomics", upgrade="never", force=T, clean=T)
```
## Guide / Overview

### Data organization and preprocessing

The library works on alignments in the form of dataframes or matrices with sites as rows and species/strains as columns. Examples:

```R
library(PhyloEpiGenomics)
data(PhyloEpiGenomics_example_data)
#nucleotide data (-> nominal scale)
head(nucl_aln)
 human chimp gorilla orangutan
1     C     C       C         C
2     T     T       T         T
3     A     A       A         A
4     G     G       G         A
5     C     C       C         C
6     T     T       T         T

#methylation data (-> interval scale)
head(meth_fraction_aln)
      human     chimp    gorilla orangutan
1 0.8500000 0.9090909 0.83333333 0.7272727
2 0.9565217 0.9565217 0.83333333 1.0000000
3 0.9500000 0.9354839 0.50000000 1.0000000
4 0.8636364 1.0000000 0.90476190 0.9259259
5 0.0000000 0.0000000 0.07142857 0.0000000
6 0.0000000 0.0000000 0.00000000 0.0000000
```

Since the library was primarily designed to analyze discretized epigenomic data, nominal data, such as nucleotides, must first be converted to integers:

```R
nucl_states_aln=sapply(nucl_aln,function(x) as.numeric(factor(x,levels=c("A","C","G","T"))))
head(nucl_states_aln)
     human chimp gorilla orangutan
[1,]     2     2       2         2
[2,]     4     4       4         4
[3,]     1     1       1         1
[4,]     3     3       3         1
[5,]     2     2       2         2
[6,]     4     4       4         4
```

Note that while this works for most implemented evolutionary models and tree reconstruction methods, it makes no sense for the "noJump" model and the adaptation of the parsimony algorithm. Therefore, the respective settings and functions should not be used if the underlying data is nominal.

If the underlying data is interval scaled like the methylation example above, it can used directly for the parsimony functions. For the other tree reconstruction functions, the interval scaled data should be discretized. The example below discretizes the [0,1] interval equally in 5 states:
```R
discretization=list(c(-0.01,0.2),c(0.2,0.4),c(0.4,0.6),c(0.6,0.8),c(0.8,1))
meth_states_aln=discretize(meth_fraction_aln,discretization)
head(meth_states_aln)
  human chimp gorilla orangutan
1     5     5       5         4
2     5     5       5         5
3     5     5       3         5
4     5     5       5         5
5     1     1       1         1
6     1     1       1         1
```
### Tree reconstruction
#### Maximum likelihood
##### Creation of evolutionary models
Tree reconstruction via maximum likelihood requires an evolutionary model. Several well-known models specific for nucleotide data are implemented: JC69 (Jukes and Cantor 1969), K80 (Kimura 1980), F81 (Felsenstein 1981), HKY85 (Hasegawa et al. 1985). For a description of those models see https://en.wikipedia.org/wiki/Models_of_DNA_evolution. In addition, our COOC (cooccurrence) model works both if the underlying data was originally nominal scaled (e.g. nucleotides or amino acids) or interval/ordinal scaled (e.g. methylation fractions). The noJump model should only be used if the underlying data was interval/ordinal scaled. See <Link to paper or bioRXiv> for the specification of the latter two models. Usage examples:

```R
#model parameters are estimated from data
my_JC69_model=make_evolutionary_model(nucl_states_aln,model="JC69")

#evolutionary models consist of a transition rate matrix Q and an equilibrium frequency vector pi
my_JC69_model
$Q
             A            C            G            T
A -0.010000000  0.003333333  0.003333333  0.003333333
C  0.003333333 -0.010000000  0.003333333  0.003333333
G  0.003333333  0.003333333 -0.010000000  0.003333333
T  0.003333333  0.003333333  0.003333333 -0.010000000

$pi
[1] 0.2975 0.1855 0.2225 0.2945

# models can also be parameterized manually
my_HKY85_model=make_evolutionary_model(model="HKY85",pi=rep(0.25,4),kappa = 2)

# make a 5 states evolutionary model for methylation data (according to the discretization used above)
my_noJump_model=make_evolutionary_model(meth_states_aln,model="noJump",nstates=length(discretization))
my_noJump_model
$Q
            [,1]        [,2]         [,3]         [,4]         [,5]
[1,] -0.00190183  0.00190183  0.000000000  0.000000000  0.000000000
[2,]  0.01165905 -0.01434051  0.002681462  0.000000000  0.000000000
[3,]  0.00000000  0.00190183 -0.010099781  0.008197951  0.000000000
[4,]  0.00000000  0.00000000  0.002681462 -0.025491611  0.022810149
[5,]  0.00000000  0.00000000  0.000000000  0.008197951 -0.008197951

$pi
[1] 0.24675 0.04025 0.05675 0.17350 0.48275
```
##### Application of evolutionary models for tree reconstruction


```R
```

```R
```