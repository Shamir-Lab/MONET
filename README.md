# MONET
 MONET is a multi-omic clustering algorithm. Recent research shows that different omics have different structure - samples that are similar in one omic are not necessarily similar in another. MONET addresses this issue using a unique formulation for multi-omic clustering. In this formulation, samples are clustered into modules, but each module is also identified by a subset of the omics. Samples in the same module are similar to one another in the omics that the module "covers". For more details of the method, see its publication: "MONET: Multi-omic module detection by omic selection".

This repository contains two parts. The first part is the code of the MONET package, which is currently available only in python. R support and more functionality are currently under development.
The second part of the repository is the code used to produce the results for MONET's paper. This code is written in R. The required libraries, data, and configuration are documented inside the code.
The repository also contains a toy multi-omic dataset to demonstrate how to use MONET.

Now follows a brief code example showing how to use MONET:
TODO.

