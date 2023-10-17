#Import Packages
#Needed for moving to output
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as pl
import itertools
from collections import Counter
import pickle

from vivarium.library.units import units, remove_units
from tumor_tcell.library.phylogeny import get_phylogeny
from tumor_tcell.experiments.main import plots_suite
from tumor_tcell.experiments.main import make_snapshot_video

#Analysis tumor-tcell modules needed
from tumor_tcell.library.data_process import data_to_dataframes
from tumor_tcell.library.data_process import control_data_to_dataframes
from tumor_tcell.library.population_analysis import division_analysis
from tumor_tcell.library.population_plots import population_plot
from tumor_tcell.library.population_plots import division_plot
from tumor_tcell.library.population_plots import death_plot

def individual_analysis(analysis_dir, experiment_id, bounds, tcells=True, lymph_nodes=False):
    # Read in the data from parent directory
    experiment_dir = analysis_dir + experiment_id
    os.chdir(experiment_dir)
    figures_out_dir = experiment_dir + '/figures'
    os.makedirs(figures_out_dir, exist_ok=True)

    #read in the data
    file_to_read = open("data_export.pkl", "rb")
    data = pickle.load(file_to_read)

    # # Check if the file exists
    # if os.path.exists("config_export.pkl"):
    #     # Write sim_config description to txt file
    #     config_to_read = open("config_export.pkl", "rb")
    #     sim_config = pickle.load(config_to_read)
    #     config_file = open("sim_config.txt", "wt")
    #     n = config_file.write(sim_config['description'])
    #     config_file.close()
    #     print(sim_config)

    #Plot the data using tumor-tcell experiment notebook and save in current directory
    fig1, fig2, fig3, fig4 = plots_suite(data, out_dir=figures_out_dir, bounds=[b * units.um for b in bounds])
    make_snapshot_video(data, bounds=[b * units.um for b in bounds], n_steps=100, out_dir=figures_out_dir)

    if tcells:
        if lymph_nodes:
            df_tumor_death, df_tcell_death, tumor_plot, tcell_plot, df_dendritic_death, dendritic_plot = data_to_dataframes(data, lymph_nodes=lymph_nodes)
        else:
            df_tumor_death, df_tcell_death, tumor_plot, tcell_plot = data_to_dataframes(data)

        # Population analysis
        divide_time_T = division_analysis(tcell_plot)
        divide_time_tumor = division_analysis(tumor_plot)

        if not divide_time_T.empty:
            division_plot(divide_data=divide_time_T, out_dir=figures_out_dir, save_name='Tcells')

        if not divide_time_tumor.empty:
            division_plot(divide_data=divide_time_tumor, out_dir=figures_out_dir, save_name='Tumors')

        population_plot(population_data=tumor_plot, cell_states=['PDL1n', 'PDL1p'], out_dir=figures_out_dir,save_name='Tumors')
        population_plot(population_data=tcell_plot, cell_states=['PD1n', 'PD1p'], out_dir=figures_out_dir,save_name='Tcells')

        if lymph_nodes:
            population_plot(population_data=dendritic_plot, cell_states=['inactive', 'active'], out_dir=figures_out_dir,
                        save_name='Dendritic')
            if not df_dendritic_death.empty:
                death_plot(death_data=df_dendritic_death, out_dir=figures_out_dir, save_name='Dendritic')

        if not df_tumor_death.empty:
            death_plot(death_data=df_tumor_death, out_dir=figures_out_dir, save_name='Tumors')
        if not df_tcell_death.empty:
            death_plot(death_data=df_tcell_death, out_dir=figures_out_dir, save_name='Tcells')


        #Export the dataframes for analysis comparison to other experiments and not need to process again
        df_tumor_death['experiment_id'] = experiment_id
        df_tcell_death['experiment_id'] = experiment_id
        tumor_plot['experiment_id'] = experiment_id
        tcell_plot['experiment_id'] = experiment_id

        if lymph_nodes:
            df_dendritic_death['experiment_id'] = experiment_id
            dendritic_plot['experiment_id'] = experiment_id

        df_tumor_death.to_csv('tumor_death.csv')
        df_tcell_death.to_csv('tcell_death.csv')
        tumor_plot.to_csv('tumor_plot.csv')
        tcell_plot.to_csv('tcell_plot.csv')

        if lymph_nodes:
            df_dendritic_death.to_csv('dendritic_death.csv')
            dendritic_plot.to_csv('dendritic_plot.csv')
    else:
        # Population analysis
        df_tumor_death, tumor_plot = control_data_to_dataframes(data)
        divide_time_tumor = division_analysis(tumor_plot)
        division_plot(divide_data=divide_time_tumor, out_dir=figures_out_dir, save_name='Tumors')
        population_plot(population_data=tumor_plot, cell_states=['PDL1n', 'PDL1p'], out_dir=figures_out_dir,
                        save_name='Tumors')
        death_plot(death_data=df_tumor_death, out_dir=figures_out_dir, save_name='Tumors')

        # Export the dataframes for analysis comparison to other experiments and not need to process again
        df_tumor_death['experiment_id'] = experiment_id
        tumor_plot['experiment_id'] = experiment_id

        df_tumor_death.to_csv('tumor_death.csv')
        tumor_plot.to_csv('tumor_plot.csv')