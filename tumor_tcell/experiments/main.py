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
    plot_agents_multigen,
    make_agent_ids,
)
from vivarium.library.units import units

# tumor-tcell imports
from tumor_tcell.experiments.control import control
# tumor-tcell composites
from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment



# simulation #1
def simulation_1(out_dir='out'):

    # experiment parameters
    total_time = 100 #0
    bounds = [1 * units.mm, 100 * units.mm]  # TODO - add units to bounds
    time_step = 60
    tumor_id = 'tumor'
    tcell_id = 'tcell'

    # get the cell configuration
    cell_config = [
        {
            'number': 3,
            'name': tumor_id,
            'type': TumorAgent,
            'config': {},
        },
        {
            'number': 1,
            'name': 'big_tumor',
            'type': TumorAgent,
            'config': {
                'tumor': {
                    'diameter': 50 * units.um,
                }
            },
        },
        {
            'number': 3,
            'name': tcell_id,
            'type': TCellAgent,
            'config': {},
        }
    ]
    make_agent_ids(cell_config)

    # configure the environment
    environment_config = {
        'type': TumorMicroEnvironment,
        'config': {
            'neighbors': {
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



# all of the experiments go here for easy access by control class
experiments_library = {
    '1': simulation_1,
}


if __name__ == '__main__':
    control(experiments_library)
