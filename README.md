# Term Paper Repository - *Dynamic Programming, Spring 2022*
This is the working repository for writing our term paper in the course [Dynamic Programming](https://kurser.ku.dk/course/a%C3%98kk08207u/) at the University of Copenhagen.

The group consists of Frederik Degn Pedersen, Christian Brauer Johanssen and Signe Holst Madsen.
## Replicating the results
To replicate the results of the paper, the interested reader would have to review the two following notebooks representing the models in A.1 and section 2 in our paper:
- [One-Asset Model](one_asset/OneAssetModel.ipynb)
- [Two-Asset Model](two_asset/TwoAssetModel.ipynb)

Each model-specific notebook is placed in a separate folder along with the modules needed for solving the model. The notebooks each contain the code for replicating the figures and computation times. 

**Figures:**

Each notebook will plot the figures associated with the given model

**Computing time**

By running the notebook the total and period-specific computation times are reported when loading the models for each solution method. 

**Model specification: Table 1 and Table 6:**

Can be found in the model initialization in [Two-Asset Model](two_asset/TwoAssetModel.py).

**Estimation:**

Fuck me
## Dependencies

We thank Jeppe Druedahl for his contributions to the consumption-saving literature and for making the code open source. Our main model draws (heavily?) from the models put forward in [this repository](https://github.com/NumEconCopenhagen/ConsumptionSavingNotebooks). 

Packages required for running the notebooks are:
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [numba](https://pypi.org/project/numba/)
- [quantecon](https://pypi.org/project/quantecon/)
- [ConSav](https://pypi.org/project/ConSav/)
- [matplotlib](https://pypi.org/project/matplotlib/)
