## Description
Gene regulatory networks play a crucial role in controlling an organism’s 
biological processes, which is why there is significant interest in developing 
computational methods that are able to extract their structure from high-throughput 
genetic data. We propose a novel efficient Bayesian method (BFCS) for discovering 
local causal relationships among triplets of (normally distributed) variables. In
our approach, we score covariance structures for each triplet in one go and incorporate
available background knowledge in the form of priors to derive posterior probabilities 
over local causal structures. Our method is flexible in the sense that it allows 
for different types of causal structures and assumptions. The proposed algorithm 
produces stable and conservative posterior probability estimates over local causal 
structures that can be used to derive an honest ranking of the most meaningful 
regulatory relationships.

## Content

The data set contains source code implementing the BFCS algorithm, which is 
described in the article titled "[A Bayesian Approach for Inferring Local Causal 
Structure in Gene Regulatory Networks](http://proceedings.mlr.press/v72/bucur18a.html)" 
by Ioan Gabriel Bucur, Tom van Bussel, Tom Claassen and Tom Heskes, as well as in
the follow-up extension titled "[Large-scale Local Causal Inference of Gene Regulatory
Relationships](https://doi.org/10.1016/j.ijar.2019.08.012). The data set also 
contains simulated data necessary for reproducing the figures in the article as 
well as routines necessary for recreating it. This research is presented in 
Chapter 3 of the PhD thesis titled "Being Bayesian about Causal Inference" by
Ioan Gabriel Bucur. The code is written in the R and C++ programming languages.

## Structure

The code is structured on the skeleton of an [R package](https://r-pkgs.org/index.html) 
package as follows:

- The folder `data` contains pre-saved simulated data, which we use to recreate
the figures from the article. The simulated data can also be reproduced using 
the `scripts/reproduce-data.R` R script. The simulated data sets are described in `R/data.R`.

- The folder `R` contains the R files that implement the BFCS algorithm and the 
routines necessary for reproducing the figures from the article. The main method
is implemented in `R/compute_BFCS_vectorized.R`.

- The folder `man` contains the documentation for the implemented functions.

- The folder `src` contains an efficient Rcpp implementation for computing the
BGe score described in [Geiger and Heckerman (2002)](https://projecteuclid.org/euclid.aos/1035844981).

- The folder `scripts` contains the script `pgm2018-article-figures.R`, which
can be run from R to produce the figures in the PGM 2018 article and the script
`ijar-article-figures.R`, which can produce the figures in the IJAR extension.
Finally `reproduce-data.R` is a script for reproducing the simulated data and
the pre-computed posterior probabilities.

- The folder `tests` contains a few basic unit tests for the R package.

- The top folder also contains the following files:
  - `DESCRIPTION` is the file describing the R package
  - `NAMESPACE` is the file specifying the functions provided by the R package
  - `LICENSE.md` is the file containing the GPL-3 license
  
## Prerequisites

In order to install the software, [R](https://cran.r-project.org/) must be 
downloaded and installed. For reproducing the tables from the articles, the TeX
system including XeTeX is required. A portable distribution satisfying this
requirement is [TeX Live](https://www.tug.org/texlive/).
  
## Installation Instructions

Download the software from GitHub with the following command:
`git clone https://github.com/igbucur/BFCS.git`. For installing and running the 
BFCS R package, several R package are required. These are specified in the 
package `DESCRIPTION` file.

To install the package, open an R instance and run (from the BFCS folder):
```
install.packages('devtools') # required package
devtools::install_deps(".", dependencies = TRUE) # install BFCS package dependencies
install.packages(".", repos = NULL, type = 'source') # install BFCS

library(BFCS) # load the package
help(package = "BFCS") # see available functions
```

## Licensing

BFCS algorithm - Bayes Factors of Covariance Structures

Copyright (C) 2020 Ioan Gabriel Bucur <ioan.gabriel.bucur@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.