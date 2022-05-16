#############
## READ ME ##
#############

The code is written and tested for 64-bit Windows machines.
It only requires free open source languages and software:


############################
## PROGRAMS AND COMPILERS ##
############################

The suggested programs and compilers are:

- ANACONDA Python 2.7 (http://docs.continuum.io/anaconda/)
- GCC Compiler with OpenMP from TDM-GCC (http://tdm-gcc.tdragon.net/download)


######################
## FOLDER STRUCTURE ##
######################

All Python files are in the main folder. 
The subfolder \cfuncs include some c-files.
The remaining folder structure is created by the program itself.


##########
## MAIN ##
##########

The model is solved and simulated by running main.py (see progress in log.txt).

It calls the following Python files in the following order:

 -> settings.py
 -> grid.py (uses check_bounds_func.c)
 -> solve.py (uses vfi.c)    
    -> vfi.c uses:
       generic_funcs.c
       base_funcs.c	
       ucon_c_func.c
       grid_search.c
       value_func.c
 -> plot_sol.py
 -> simulate.py (uses sim_func.c)
    -> sim_func.c uses:
       generic_funcs.c
       base_funcs.c
 -> plot_sim.py
 -> plot_welfare.py
 -> plot_robustnes.py

Explanations of each file is included in its header.


#############################
## Standalone Python-files ##
#############################

 -> nonconvexchoiceset.py (produces figures of the choice set)