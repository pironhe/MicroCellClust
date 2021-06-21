# MicroCellClust

This repository contains a new version of the implementation of the MicroCellClust optimisation problem. It is designed to search for rare cell-clusters (typically less than 10% of the number of cells) that present a highly specific gene expression in scRNA-seq data. It has been implemented in the Scala programming language. This repository also contains an interface to run MicroCellClust in R. This new version tunes the meta-parameters automatically.

Version 1.2 is also available on Zenodo http://doi.org/10.5281/zenodo.4580332



## Using MicroCellClust in R

### RScala

To run MicroCellClust in R, the user must install the `rscala` package (https://cran.r-project.org/web/packages/rscala/index.html), that can call the Scala implementation of MicroCellClust from within R. Once installed, he can check wether there is a version of Scala on his computer using:

``` R
rscala::scalaConfig()
```
This function tries to find Scala and Java on the user's computer and, if needed, downloads and installs Scala and Java in the user's `~/.rscala` directory.


### Initializing MicroCellClust

To run MicroCellClust, the user must first set his working directory (`setwd`) to the location of the `MicroCellClust`directory, and then load the functions contained in the `micro_cell_clust.R` script.
``` R
source("./src/main/R/micro_cell_clust.R")
```
Then, he must initialize an instance of RScala containing the MicroCellClust source code.
``` R
mcc.rs = initMccRscala()
```
This `mcc.rs` object needs to be passed as the first argument to any R function using Scala code.


### Running MicroCellClust with automatic tuning

Running MicroCellClust with automatic tuning is done using the following function:

``` R
results = runAutomaticTuning(mcc.rs, data)
```

Since the parameters are automatically tuned, it is not necessary to pass them as argument. It is still possible to pass values of `kappa` and `mu` by using this function:

``` R
results = runMCC(mcc.rs, data, kappa = 1, nNeg= 0.1)
```

But the expression matrix requires to be preprocessed manually, while it is done automatically in the `runAutomaticTuning` function.

###Example

An example can be found in the file called `usage_example.R` at the root of this repository
