### Regularized NCA for estimating TF activity

[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![Generic badge](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/Roy-lab/EstimateNCA/releases/tag/v1.0.0)

Regularized NCA uses the same framework as original NCA ([Liao et al. PNAS 2003](https://doi.org/10.1073/pnas.2136632100)), but uses a modified LASSO formulation to incorporate edge confidence from input network.
Briefly, NCA uses a two step iterative method, estimates TFA profiles from current network, and then estimate regression coefficients of network from the new TFA, and repeats until convergance.
In regularized TFA, when estimating regression coefficients, we use a modified LASSO to incorporate edge confidence in order to shrink low confidence interactions from the model.

We add the resulting TFA profiles to our inference method [MERLIN-P](https://github.com/Roy-lab/merlin-p) ([Roy et al. PLOS Comput Biol 2013](https://doi.org/10.1371/journal.pcbi.1003252), [Siahpirani & Roy NAR 2017](https://doi.org/10.1093/nar/gkw963)) which we call MERLIN-P+TFA.

![alt text](example/tfa_overview.png "Overview of MERLIN-P+TFA. We start with an expression matrix and an input prior network. TF activity profile is estimated using regularized NCA, and final inferred network in inferred using estimated TFA and the input expression matrix and the prior network.")

### How to use

In order to use the program, you will need:

* a prior network with edge confidence (first column is TF, second column is target gene, third column is confidence, see example under example/in). The edge weights could be between 0 and 1, lower (near 0) means low confidence.
   * This could be from any source (ChIP, motif, etc.). For mammalian cell lines we have used [PIQ](http://piq.csail.mit.edu/) that can integrate DNase/ATAC with a motif network.
* a expression matrix (first column is gene names, the rest are expression values, no header, see example under example/in).

You will need to estimate TFA multiple time (rand inits) and take average (see scripts under example/stability_example/). For yeast dataset we had used 10 rand inits, for mammalian data, we used 100 rand inits.

In order to use the estimated TFA profile in network inference process, we add a suffix to TF name so we would be able to concatenate the TFA profile to expression matrix and use bot expression and TFA profile of TFs.
```
paste <(cut -f1 tfa.txt |awk '{printf("%s_nca\n",$1)}' ) <(cut -f2- tfa.txt ) > tfa_with_suffix.txt
```

### Run

You will need GSL library to compile the code. In order to compile to code, navigate to code directory and make.

You can run the program as:

```
./NCALearner -d example/in/exp.txt -r example/in/tfs.txt -g example/in/targets.txt -p example/in/prior.txt -o example/out/
```

The input options are: 
```
-h  prints the usage.
-d	is the expression matrix.
-r  is the list of TFs (we only consider these).
-g  is the list of targets (we only consider these).
-p  is the prior network.
-l  is the lambda for NCA.
-o  is output directory.
```

We have some example input/output files in:

```
example/in:
	exp.txt
	prior.txt
	targets.txt
	tfs.txt

example/out:
	adj.txt
	tfa.txt
```

You can find examples on how to take average of multiple rand inits and compare them to each other under example/stability_example/ 

![alt text](example/stability_example/corr_stability.png "Absolute value of correlation of estimated TF activities between different rand inits. Each profile is average over 10 rand inits. 10 such profiles were generated and were compared to each other. We show 25%, 50%, and 75% quantiles of the correlation values.")
