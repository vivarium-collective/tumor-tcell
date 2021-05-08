import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl

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
        pl.savefig(out_dir+'/'+save_name+'_division.png', format='png', bbox_inches='tight', dpi=300)

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

    pl.figure(figsize=(8, 4))
    ttl = sns.lineplot(data=total_cell, x="time", y='cell', label='total')
    ttl_state1 = sns.lineplot(data=state_1, x="time", y='cell', label=cell_states[0])
    ttl_state2 = sns.lineplot(data=state_2, x="time", y='cell', label=cell_states[1])

    pl.title("Total "+save_name)
    pl.legend(title="Cell type")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_total.png', format='png', bbox_inches='tight', dpi=300)


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
    death_cell = sns.lineplot(data=death_plot, x="time", y='death count', hue='death type')
    pl.title("# of deaths")
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if save_name is not None:
        pl.savefig(out_dir + '/' + save_name + '_death.png', format='png', bbox_inches='tight', dpi=300)