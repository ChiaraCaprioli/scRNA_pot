# TO DO:
# Fix size in MCMC diagnostic plots
# Add table with acceptance rate by test/parameters

## Setup
import importlib
import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import arviz as az
import tensorflow as tf

## Run scCODA
def RunScCODA(path_cell_counts, path_save, contrasts, group_var, alpha, seed):
    
    ## Set path
    PATH_scCODA = os.path.join(path_save, "scCODA")
    os.makedirs(os.path.join(path_save, "scCODA"), exist_ok = True) 

    ## Set seed
    tf.random.set_seed(seed = seed)

    ## Prepare data
    # Load data
    cell_counts = pd.read_csv(path_cell_counts, index_col=0)
    cell_counts = cell_counts.drop(['total_cell_count'], axis=1)

    # Convert to anndata
    data_all = dat.from_pandas(cell_counts, covariate_columns=["sample_id", group_var])

    ## Compositional analysis by contrast
    for i in contrasts:
        
        # Set path
        PATH_SAVE = os.path.join(PATH_scCODA, i)
        os.makedirs(os.path.join(PATH_scCODA, i), exist_ok = True) 

        # Subset anndata to contrast of interest
        data_contrast = data_all[data_all.obs[group_var].isin([i.split('_')[0], i.split('_')[1]])]

        # Set reference
        ref = i.split('_')[1]

        # Set model
        model = mod.CompositionalAnalysis(
            data_contrast, formula=f"C({group_var}, Treatment('{ref}'))", 
            reference_cell_type = "automatic"
            #automatic_reference_absence_threshold = 0.05 # default = 0.05
            ) 
            
        # Run MCMC
        res = model.sample_hmc() # desired acceptance rate 0.4-0.9 (otherwise check sampling issues)

        # Set desired FDR
        res.set_fdr(est_fdr = alpha)

        # Plot and save MCM diagnostics
        az.plot_trace(
            res,
            divergences=False,
            var_names=["alpha", "beta"],
            coords={"cell_type": res.posterior.coords["cell_type_nb"]},
        )
        plt.savefig(os.path.join(PATH_SAVE, f"{i}_MCMC_diagn.png"))
        
        # Save all results as pickle 
        ## saving
        PATH_PICKLE = os.path.join(PATH_SAVE, f"{i}.pkl")
        res.save(PATH_PICKLE)

        ## loading
        with open(PATH_PICKLE, "rb") as f:
            res = pkl.load(f)

        # Save summary as csv
        res.effect_df.to_csv(os.path.join(PATH_SAVE, f"{i}.csv"))
        
