# MONET
 MONET is a multi-omic clustering algorithm. Recent research shows that different omics have different structure - samples that are similar in one omic are not necessarily similar in another. MONET addresses this issue using a unique formulation for multi-omic clustering. In this formulation, samples are clustered into modules, but each module is also identified by a subset of the omics. Samples in the same module are similar to one another in the omics that the module "covers". For more details of the method, see its publication: "MONET: Multi-omic module detection by omic selection".

This repository contains two parts. The first part is the code of the MONET package, which is currently available only in python. R support and more functionality are currently under development. It is places under the directory monet/.
The second part of the repository is the code used to produce the results for MONET's paper. This code is written in R. The required libraries, data, and configuration are documented inside the code. It is places under the directory R_code/.
The repository also contains a toy multi-omic dataset to demonstrate how to use MONET. It is places under the directory data/.

Now follows a brief code example showing how to use MONET.
```{python}
from monet import data_prep
from monet.monet import Monet

# helper function to read the example csv files.
path_to_csv_dir = '/path/to/data/'
sim_mats, _ = data_prep.process_data(path=path_to_csv_dir)
# main_loop is the main Monet function. See its documentation for a full explanation on the parameters and return values.
monet_ret = Monet().main_loop(data=sim_mats, num_of_samples_in_seed=26, num_of_seeds=15, min_mod_size=13)
monet_results, total_time, super_g, iteration_times, total_weight = monet_ret
#all_modules is now a dict mapping module names to module objects.
all_modules = monet_results.modules
arbitrary_module = list(all_modules.values())[0]
# arbitrary_module now contains a dictionary of its patients (from id to Patient object)
# and a dictionary of its omics (from id to Omic object)
arbitrary_module.patients
arbitrary_module.omics
```

