# Term Paper Repository - *Dynamic Programming, Spring 2022*
This is the working repository for writing our term paper in the course [Dynamic Programming](https://kurser.ku.dk/course/a%C3%98kk08207u/) at the University of Copenhagen.

The group consists of Frederik Degn Pedersen, Christian Brauer Johanssen and Signe Holst Madsen.

## Replicating the results
To replicate the results of the paper, the interested reader would have to review the two following notebooks representing the models in A.1 and section 2 in our paper:
- [One-Asset Model](one_asset/OneAssetModel.ipynb)
- [Two-Asset Model](two_asset/TwoAssetModel.ipynb)

Each model-specific notebook is placed in a separate folder along with the modules needed for solving the model. The notebooks each contain the code for replicating the figures and computation times. 

**Model specification: Table 1 and Table 6:**

Can be found in the model initialization in [Two-Asset Model module](two_asset/TwoAssetModel.py). 

**Computing time**

By running the notebook, the total and period-specific computation times are reported when loading the models for each solution method. These timings make up Table 2 in our paper. 

**MPCs:**

The average MPCs in Table 3 and MPCs from Table 4 sensitivity analysis are also present in notebook. As default, the MPCs are cross-computed in the simulation part, but to create the non cross-computed MPCs a boolean for ``cross_compute=False`` can be set when initiating the model.

**Figures:**

Each notebook will plot the figures associated with the given model.

**Estimation:**
To illustrate that we are indeed able to estimate the "true" parameters. Nothing fancy going on there.. 
- [Two-Asset Estimation Notebook](two_asset/estimation.ipynb)

## Dependencies

We thank Jeppe Druedahl for his contributions to the consumption-saving literature and for making the code open source. Our main model draws (heavily?) from the models put forward in [this repository](https://github.com/NumEconCopenhagen/ConsumptionSavingNotebooks). 

Packages required for running the notebooks are:
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [numba](https://pypi.org/project/numba/)
- [quantecon](https://pypi.org/project/quantecon/)
- [ConSav](https://pypi.org/project/ConSav/)
- [matplotlib](https://pypi.org/project/matplotlib/)
