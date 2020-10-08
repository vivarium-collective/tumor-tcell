"""
========================
Tumor/T-cell Experiments
========================

Experiments can be triggered from the command line:

```
$ python tumor_tcell/experiments/main.py [experiment_name]
```
"""

# vivarium-core imports
from vivarium.core.composition import (
    agent_environment_experiment,
    make_agent_ids,
    EXPERIMENT_OUT_DIR,
)
from vivarium.library.units import units, remove_units
from vivarium.core.control import Control

# plots
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_cell.plots.multibody_physics import plot_tags, plot_snapshots

# tumor-tcell imports
from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment


out_dir = EXPERIMENT_OUT_DIR
TIMESTEP = 60


# simulation #1
def simulation_1():

    # experiment parameters
    total_time = 1000
    bounds = [200 * units.um, 200 * units.um]
    time_step = 60
    tumor_id = 'tumor'
    tcell_id = 'tcell'

    # get the cell configuration
    cell_config = [
        {
            'number': 1,
            'name': tumor_id,
            'type': TumorAgent,
            'config': {
                'time_step': TIMESTEP,
            },
        },
        {
            'number': 1,
            'name': tcell_id,
            'type': TCellAgent,
            'config': {
                'time_step': TIMESTEP,
            },
        }
    ]
    make_agent_ids(cell_config)

    # configure the environment
    environment_config = {
        'type': TumorMicroEnvironment,
        'config': {
            'neighbors_multibody': {
                'time_step': TIMESTEP,
                'bounds': bounds,
            }
        },
    }

    # make the experiment using helper function agent_environment_experiment
    experiment = agent_environment_experiment(
        agents_config=cell_config,
        environment_config=environment_config
    )

    # run the simulation
    experiment.update(total_time)
    experiment.end()

    # retrieve data from the emitter
    data = experiment.emitter.get_data()

    # separate out tcell and tumor data for multigen plots
    tcell_data = {}
    tumor_data = {}
    for time, time_data in data.items():
        all_agents_data = time_data['agents']
        tcell_data[time] = {
            'agents': {
                agent_id: agent_data
                for agent_id, agent_data in all_agents_data.items()
                if tcell_id in agent_id}}
        tumor_data[time] = {
            'agents': {
                agent_id: agent_data
                for agent_id, agent_data in all_agents_data.items()
                if tumor_id in agent_id}}

    # make multigen plot for tcells and tumors
    plot_settings = {}
    plot_agents_multigen(tcell_data, plot_settings, out_dir, tcell_id)
    plot_agents_multigen(tumor_data, plot_settings, out_dir, tumor_id)


    # snapshots plot
    # extract data
    multibody_config = remove_units(environment_config['config']['neighbors_multibody'])
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    # fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_data = {
        'agents': agents,
        'fields': {},
        'config': multibody_config}
    plot_config = {
        'fields': [],
        'n_snapshots': 8,
        'agent_shape': 'circle',
        'out_dir': out_dir}
    plot_snapshots(plot_data, plot_config)

    print()
    print('saved at: {}'.format(out_dir))



# all of the experiments go here for easy access by control class
experiments_library = {
    '1': simulation_1,
}


if __name__ == '__main__':
    Control(
        experiments=experiments_library
    )