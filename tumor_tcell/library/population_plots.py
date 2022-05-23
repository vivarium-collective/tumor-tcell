import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
import os
import numpy as np

def division_plot(divide_data, out_dir = None, save_name = None):

    #Import Plot settings
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


    #Plot number of T cell divisions
    pl.figure(figsize=(8, 4))
    div_cell_T = sns.lineplot(data=divide_data, x="time", y='total_division')
    pl.title("# of divisions")
    if save_name is not None:
        pl.savefig(out_dir+'/'+save_name+'_division.png', transparent=True, format='png', bbox_inches='tight', dpi=300)

def population_plot(population_data, cell_states, out_dir=None, save_name=None):
    # Import Plot settings
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Plot total cells and see how changing over time
    total_cell = population_data.groupby('time')['cell'].nunique().reset_index()
    cell_state_df = population_data.groupby(['time', 'cell_state'])['cell'].nunique().reset_index()
    state_1 = cell_state_df.loc[cell_state_df['cell_state'] == cell_states[0]]
    state_2 = cell_state_df.loc[cell_state_df['cell_state'] == cell_states[1]]

    # reset index
    total_cell.reset_index(inplace=True, drop=True)
    state_1.reset_index(inplace=True, drop=True)
    state_2.reset_index(inplace=True, drop=True)

    pl.figure(figsize=(8, 4))
    ttl = sns.lineplot(data=total_cell, x="time", y='cell', label='total')
    ttl_state1 = sns.lineplot(data=state_1, x="time", y='cell', label=cell_states[0])
    ttl_state2 = sns.lineplot(data=state_2, x="time", y='cell', label=cell_states[1])

    pl.title("Total "+save_name)
    pl.legend(title="Cell type")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_total.png', transparent=True, format='png', bbox_inches='tight', dpi=300)


def death_plot(death_data, out_dir=None, save_name=None):
    # Import Plot settings
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    #Plot total number of deaths and type
    total_col = [col for col in death_data.columns if 'total' in col]
    pl.figure(figsize=(8, 4))
    death_plot = pd.melt(death_data, id_vars= ['death', 'time'], value_vars= total_col)
    death_plot.rename(columns={'variable':'death type', 'value' : 'death count'}, inplace=True)

    # reset index
    death_plot.reset_index(inplace=True, drop=True)

    death_cell = sns.lineplot(data=death_plot, x="time", y='death count', hue='death type')
    pl.title("# of deaths")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_death.png', transparent=True, format='png', bbox_inches='tight', dpi=300)


def death_group_plot(death_plot_list, out_dir=None, save_name=None):
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    death_plot_l = []
    for experiment in death_plot_list:
        # Plot total number of deaths and type
        total_col = [col for col in experiment.columns if 'total' in col]
        if len(experiment) > 1:
            death_plot = pd.melt(experiment, id_vars=['death', 'time', 'experiment_name'], value_vars=total_col)
            death_plot.rename(columns={'variable': 'death type', 'value': 'death count'}, inplace=True)
            death_plot_l.append(death_plot)

    if not death_plot_l:
        raise Warning('no death in these experiments, skipping death_group_plot')
        return

    # Concatenate all
    death_plot = pd.concat(death_plot_l)

    # Separate out total death
    total_death_plot = death_plot.loc[death_plot['death type'] == 'total_death']
    other_death_plot = death_plot.loc[~(death_plot['death type'] == 'total_death')]

    # reset index
    total_death_plot.reset_index(inplace=True, drop=True)
    other_death_plot.reset_index(inplace=True, drop=True)

    # Plot figures
    pl.figure(figsize=(8, 4))
    death_cell = sns.lineplot(data=total_death_plot, x="time", y='death count', hue='experiment_name')
    pl.title("# of deaths")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_death.png', transparent=True, format='png', bbox_inches='tight',
                   dpi=300)

    pl.figure(figsize=(8, 4))
    death_cell = sns.lineplot(data=other_death_plot, x="time", y='death count', hue='experiment_name',
                              style='death type')
    pl.title("# of deaths")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_death_subtypes.png', transparent=True, format='png',
                   bbox_inches='tight', dpi=300)


def population_group_plot(cell_plot_list, cell_states, out_dir=None, save_name=None):
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Plot total cells and see how changing over time
    experiment_plot_list = []
    cell_state_1_list = []
    cell_state_2_list = []

    for experiment in cell_plot_list:
        # Get all for total comparisons
        total_cell = experiment.groupby(['time', 'experiment_name'])['cell'].nunique().reset_index()
        experiment_plot_list.append(total_cell)

        # Get all for individual comparisons
        cell_state_df = experiment.groupby(['time', 'cell_state', 'experiment_name'])['cell'].nunique().reset_index()
        state_1 = cell_state_df.loc[cell_state_df['cell_state'] == cell_states[0]]
        cell_state_1_list.append(state_1)
        state_2 = cell_state_df.loc[cell_state_df['cell_state'] == cell_states[1]]
        cell_state_2_list.append(state_2)

    # Concatenate all
    experiment_plot = pd.concat(experiment_plot_list)
    cell_state_plot_1 = pd.concat(cell_state_1_list)
    cell_state_plot_2 = pd.concat(cell_state_2_list)
    cell_state_all = pd.concat([cell_state_plot_1, cell_state_plot_2])

    # reset index
    experiment_plot.reset_index(inplace=True, drop=True)
    cell_state_all.reset_index(inplace=True, drop=True)

    # Create plot
    pl.figure(figsize=(8, 4))
    ttl_1 = sns.lineplot(data=experiment_plot, x="time", y='cell', hue='experiment_name')
    pl.title("Total " + save_name)
    pl.legend(title="Experiment")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        save_path = os.path.join(out_dir, save_name + '_total.png')
        pl.savefig(save_path, transparent=True, format='png', bbox_inches='tight',
                   dpi=300)

    # Create plot
    pl.figure(figsize=(8, 4))
    ttl_2 = sns.lineplot(data=cell_state_all, x="time", y='cell', hue='experiment_name', style='cell_state')
    pl.title("Total Subtype of " + save_name)
    pl.legend(title="Experiment")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        save_path = os.path.join(out_dir, save_name + '_total_subtype.png')
        pl.savefig(save_path, transparent=True, format='png',
                   bbox_inches='tight', dpi=300)

def cytotoxicity_group_plot(
        cell_plot_list, exp_1, cntrl_1, exp_2, cntrl_2,
        out_dir=None, save_name=None
):
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Combine data for comparisons
    experiment_plot_list = []
    for experiment in cell_plot_list:
        total_cell = experiment.groupby(['time', 'experiment_name'])['cell'].nunique().reset_index()
        experiment_plot_list.append(total_cell)
    experiment_plot = pd.concat(experiment_plot_list)

    # Reshape the data for calculating cytotoxicity with controls
    tt = experiment_plot.pivot(index='time', columns='experiment_name', values='cell').reset_index()
    exp_1 = exp_1
    cntrl_1 = cntrl_1
    exp_2 = exp_2
    cntrl_2 = cntrl_2
    col_1 = ['time', exp_1, cntrl_1]
    col_2 = ['time', exp_2, cntrl_2]
    tt_1 = tt[col_1]
    tt_2 = tt[col_2]

    # calculate cytotoxicity with controls
    tt_1['cytotoxicity'] = (tt_1[cntrl_1] - tt_1[exp_1]) / tt_1[cntrl_1] * 100
    tt_2['cytotoxicity'] = (tt_2[cntrl_2] - tt_2[exp_2]) / tt_2[cntrl_2] * 100

    # reshape data for plotting together
    tm_1 = pd.melt(tt_1, id_vars='time', value_vars='cytotoxicity')
    tm_1['experiment_name'] = exp_1
    tm_1.rename(columns={'value': 'cytotoxicity'}, inplace=True)
    tm_2 = pd.melt(tt_2, id_vars='time', value_vars='cytotoxicity')
    tm_2['experiment_name'] = exp_2
    tm_2.rename(columns={'value': 'cytotoxicity'}, inplace=True)
    cytotoxic_plot = pd.concat([tm_1, tm_2])

    # reset index
    cytotoxic_plot.reset_index(inplace=True, drop=True)

    # Create plot
    pl.figure(figsize=(8, 4))
    ttl_1 = sns.lineplot(data=cytotoxic_plot, x="time", y='cytotoxicity', hue='experiment_name')
    pl.title("Total " + save_name)
    pl.legend(title="Experiment")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        save_path = os.path.join(out_dir, save_name + '_cytotoxicity.png')
        pl.savefig(save_path, transparent=True, format='png', bbox_inches='tight', dpi=300)

    return cytotoxic_plot


def cytotoxicity_rep_plot(
        num_rep,
        save_name='cytotoxicity_rep_plot',
        analysis_dir='out/killing_experiments/'
):
    #settings for the plotting
    SMALL_SIZE = 18
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 24

    pl.rc('font', size=SMALL_SIZE)  # controls default text sizes
    pl.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Get csv saved in experiment id library
    save_dir = analysis_dir + 'Multiple_killing_analysis/'
    analysis_out_dir = save_dir + save_name + '_n55_rep'  #####################
    os.makedirs(analysis_out_dir, exist_ok=True)

    experiment_list = []
    for n in range(1, num_rep + 1, 1):
        s = 'MHCI_Reduction_' + save_name + '_ncells_55_exp' + str(n)
        experiment_list.append(s)

    # Create a list of dataframes for concatenation
    high_cytotoxicity_list = []

    # read in dataframes
    for experiment in experiment_list:
        experiment_id = experiment
        experiment_dir = analysis_dir + experiment_id + '/killing_PDL1_50_PDL1_0/'

        df_cytotoxicity = pd.read_csv(experiment_dir + experiment + '_cytotoxicity.csv', index_col=0)
        df_cytotoxicity['replicate'] = experiment
        high_cytotoxicity_list.append(df_cytotoxicity)

    # concantenate columns and calculate mean and sem
    df = pd.concat(high_cytotoxicity_list)
    dfb = df.groupby(['time', 'experiment_name']).agg({'cytotoxicity': ['mean', 'sem']})
    dfb.columns = dfb.columns.droplevel(0)
    dfb = dfb.rename_axis(None, axis=1)
    dfb.reset_index(inplace=True)

    pl.figure(figsize=(8, 4))
    for experiment in dfb.experiment_name.unique():
        df_plot = dfb[dfb.experiment_name == experiment]

        # create vectors to use for SEM plotting
        M_new_vec = np.array(df_plot['mean'])
        Sigma_new_vec = np.array(df_plot['sem'])
        lower_bound = M_new_vec - Sigma_new_vec
        upper_bound = M_new_vec + Sigma_new_vec

        # reset index
        df_plot.reset_index(inplace=True, drop=True)

        # Create plot
        ttl_1 = sns.lineplot(data=df_plot, x="time", y='mean')
        pl.title("Total " + save_name)
        pl.ylabel('cytotoxicity')
        pl.fill_between(df_plot['time'], lower_bound, upper_bound, alpha=.3)

    #save the plot
    if save_name is not None:
        pl.savefig(analysis_out_dir + '/' + save_name + '_cytotoxicity.png', transparent=True, format='png',
                   bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    cytotoxicity_rep_plot(num_rep=2, save_name='2_rep')