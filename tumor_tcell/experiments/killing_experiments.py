#Import important packages for testing out the models throughout this jupyter notebook
#External modules needed
import seaborn as sns
import matplotlib.pyplot as pl
import pandas as pd
import os

#Experiment tumor-tcell modules needed
from tumor_tcell.experiments.main import (tumor_tcell_abm, plots_suite, get_tcells, get_tumors)
from vivarium.library.units import units, remove_units

#Analysis tumor-tcell modules needed
from tumor_tcell.library.data_process import data_to_dataframes
from tumor_tcell.library.data_process import control_data_to_dataframes

#Analysis tumor-tcell modules needed
from tumor_tcell.library.population_plots import death_group_plot
from tumor_tcell.library.population_plots import population_group_plot
from tumor_tcell.library.population_plots import cytotoxicity_group_plot


def killing_experiment(
        N_cells, save_name, num_rep=None, tumor_t_ratio=1,
        bounds=[250 * units.um, 250 * units.um],
        total_time=48600
):
    #run each replicate
    for n in range(1, num_rep + 1, 1):

        experiment_name = f'MHCI_Reduction_{save_name}_ncells_{str(N_cells)}_exp{str(n)}'
        parent_dir = 'out/killing_experiments/'

        # Make a new folder for each experiment
        out_dir = os.path.join(parent_dir, experiment_name)
        exp_out_dir_1 = os.path.join(out_dir, 'PDL1n_100percent_1')
        exp_out_dir_2 = os.path.join(out_dir, 'PDL1n_50percent_1')
        exp_out_dir_3 = os.path.join(out_dir, 'PDL1n_50percent_Cntrl_1')
        exp_out_dir_4 = os.path.join(out_dir, 'PDL1n_100percent_Cntrl_1')
        exp_out_dir_5 = os.path.join(out_dir, 'killing_PDL1_50_PDL1_0')
        for exp_dir in [
            parent_dir, out_dir, exp_out_dir_1, exp_out_dir_2, exp_out_dir_3, exp_out_dir_4, exp_out_dir_5]:
            os.makedirs(exp_dir, exist_ok=True)


        ##Experiment to compare to killing data 1:1 with 0 PDL1+ tumors
        N_TUMORS = N_cells
        N_TCELLS = int(N_cells/tumor_t_ratio)
        relative_pd1n = 0.25
        relative_pdl1n = 1
        DEFAULT_TUMORS = get_tumors(number=N_TUMORS, relative_pdl1n=relative_pdl1n)
        DEFAULT_TCELLS = get_tcells(number=N_TCELLS, relative_pd1n=relative_pd1n)

        # global parameters
        TIMESTEP = 60
        BOUNDS = bounds
        data = tumor_tcell_abm(
            total_time=total_time,
            tumors=DEFAULT_TUMORS,
            tcells=DEFAULT_TCELLS,
            halt_threshold=500,
            emit_step=60,
            bounds=BOUNDS,
            sim_step=10 * TIMESTEP,
        )
        data = remove_units(data)

        # Plot the data using tumor-tcell experiment notebook and save in current directory
        fig1, fig2, fig3 = plots_suite(data, out_dir=exp_out_dir_1, bounds=[b * units.um for b in BOUNDS])

        # Extract dataframes for plotting
        df_tumor_death_1, df_tcell_death_1, tumor_plot_1, tcell_plot_1 = data_to_dataframes(data)



        ##Experiment to compare to killing data 1:1 with 50% PDL1+ tumors
        relative_pdl1n = 0.5
        DEFAULT_TUMORS = get_tumors(number=N_TUMORS, relative_pdl1n=relative_pdl1n)
        DEFAULT_TCELLS = get_tcells(number=N_TCELLS, relative_pd1n=relative_pd1n)

        data = tumor_tcell_abm(
            total_time=total_time,
            tumors=DEFAULT_TUMORS,
            tcells=DEFAULT_TCELLS,
            halt_threshold=500,
            emit_step=60,
            bounds=BOUNDS,
            sim_step=10 * TIMESTEP,
        )
        data = remove_units(data)

        # Plot the data using tumor-tcell experiment notebook and save in current directory
        fig1, fig2, fig3 = plots_suite(data, out_dir=exp_out_dir_2, bounds=[b * units.um for b in BOUNDS])

        # Extract data for plotting
        df_tumor_death_2, df_tcell_death_2, tumor_plot_2, tcell_plot_2 = data_to_dataframes(data)



        ##Experiment to compare to killing data 0:1 with 50% PDL1+ tumors
        ##Control
        data = tumor_tcell_abm(
            bounds=BOUNDS,
            n_tumors=N_cells,
            n_tcells=0,
            tumors_state_PDL1n=0.5,
            tcells_state_PD1n=relative_pd1n,
            total_time=total_time,
            sim_step=10 * TIMESTEP,
            halt_threshold=500,
            time_step=TIMESTEP,
            emit_step=60,
        )
        data = remove_units(data)

        # Plot the data using tumor-tcell experiment notebook and save in current directory
        fig1, fig2, fig3 = plots_suite(data, out_dir=exp_out_dir_3, bounds=[b * units.um for b in BOUNDS])

        # Extract data for plotting
        df_tumor_death_3, tumor_plot_3 = control_data_to_dataframes(data)



        ##Experiment to compare to killing data 0:1 with 0% PDL1+ tumors
        ##Control
        data = tumor_tcell_abm(
            bounds=BOUNDS,
            n_tumors=N_cells,
            n_tcells=0,
            tumors_state_PDL1n=1,
            tcells_state_PD1n=relative_pd1n,
            total_time=total_time,
            sim_step=10 * TIMESTEP,
            halt_threshold=500,
            time_step=TIMESTEP,
            emit_step=60, )
        data = remove_units(data)

        # Plot the data using tumor-tcell experiment notebook and save in current directory
        fig1, fig2, fig3 = plots_suite(data, out_dir=exp_out_dir_4, bounds=[b * units.um for b in BOUNDS])

        # Extract data for plotting
        df_tumor_death_4, tumor_plot_4 = control_data_to_dataframes(data)


        ###Analyze all together for cytotoxicity
        # Add experimental tag to each dataframe
        name_exp_1 = '0% PDL1+'
        name_exp_2 = '50% PDL1+'
        name_exp_3 = '50% PDL1+ Cntrl'
        name_exp_4 = '0% PDL1+ Cntrl'

        df_tcell_death_1['experiment_name'] = name_exp_1
        df_tumor_death_1['experiment_name'] = name_exp_1
        tcell_plot_1['experiment_name'] = name_exp_1
        tumor_plot_1['experiment_name'] = name_exp_1

        df_tcell_death_2['experiment_name'] = name_exp_2
        df_tumor_death_2['experiment_name'] = name_exp_2
        tcell_plot_2['experiment_name'] = name_exp_2
        tumor_plot_2['experiment_name'] = name_exp_2

        df_tumor_death_3['experiment_name'] = name_exp_3
        tumor_plot_3['experiment_name'] = name_exp_3

        df_tumor_death_4['experiment_name'] = name_exp_4
        tumor_plot_4['experiment_name'] = name_exp_4

        # Combine together in list form for analysis together
        df_tcell_death_list = [df_tcell_death_1, df_tcell_death_2]
        df_tumor_death_list = [df_tumor_death_1, df_tumor_death_2, df_tumor_death_3, df_tumor_death_4]
        tcell_plot_list = [tcell_plot_1, tcell_plot_2]
        tumor_plot_list = [tumor_plot_1, tumor_plot_2, tumor_plot_3, tumor_plot_4]

        # save experiment dataframe in case need for later
        # os.chdir(exp_out_dir_5)

        tcell_death = pd.concat(df_tcell_death_list)
        tumor_death = pd.concat(df_tumor_death_list)
        tumor = pd.concat(tumor_plot_list)
        tcell = pd.concat(tcell_plot_list)

        tumor_death.to_csv(exp_out_dir_5 + 'tumor_death.csv')
        tcell_death.to_csv(exp_out_dir_5 + 'tcell_death.csv')
        tumor.to_csv(exp_out_dir_5 + 'tumor.csv')
        tcell.to_csv(exp_out_dir_5 + 'tcell.csv')

        cytotoxicity = cytotoxicity_group_plot(cell_plot_list=tumor_plot_list, exp_1=name_exp_1, cntrl_1=name_exp_4,
                                               exp_2=name_exp_2, cntrl_2=name_exp_3,
                                               out_dir=exp_out_dir_5, save_name=experiment_name)
        cytotoxicity.to_csv(exp_out_dir_5 + experiment_name + '_cytotoxicity.csv')

        death_group_plot(death_plot_list=df_tcell_death_list, out_dir=exp_out_dir_5, save_name='Tcells')
        death_group_plot(death_plot_list=df_tumor_death_list, out_dir=exp_out_dir_5, save_name='Tumors')

        population_group_plot(cell_plot_list=tcell_plot_list, cell_states=['PD1n', 'PD1p'], out_dir=exp_out_dir_5,
                              save_name='Tcells')
        population_group_plot(cell_plot_list=tumor_plot_list, cell_states=['PDL1n', 'PDL1p'], out_dir=exp_out_dir_5,
                              save_name='Tumor')


if __name__ == '__main__':
    killing_experiment(N_cells=20, save_name='4', num_rep=2,
                       total_time=1000
                       )
