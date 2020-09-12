"""
==============
Experiments
==============
"""

import os
import uuid
import argparse

from vivarium.core.composition import (
    agent_environment_experiment,
    plot_agents_multigen,
)

from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment

# directories
from tumor_tcell import EXPERIMENT_OUT_DIR


# helper functions
def make_agent_ids(agents_config):
    agent_ids = []
    for config in agents_config:
        number = config.get('number', 1)
        if 'name' in config:
            name = config['name']
            if number > 1:
                new_agent_ids = [name + '_' + str(num) for num in range(number)]
            else:
                new_agent_ids = [name]
        else:
            new_agent_ids = [str(uuid.uuid1()) for num in range(number)]
        config['ids'] = new_agent_ids
        agent_ids.extend(new_agent_ids)
    return agent_ids



def simulation_1(
    out_dir='out'
):
    total_time = 1000
    time_step = 60
    tumor_id = 'tumor'
    tcell_id = 'tcell'

    # configure the cells
    cell_config = [{
        'number': 3,
        'name': tumor_id,
        'type': TumorAgent,
        'config': {
            'time_step': time_step,
        },
        },
        {
        'number': 3,
        'name': tcell_id,
        'type': TCellAgent,
        'config': {
            'time_step': time_step,
        },
        }]
    make_agent_ids(cell_config)

    # configure the environment
    environment_config = {
        'type': TumorMicroEnvironment,
        'config': {
            'neighbors': {
                'time_step': time_step,
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

    # retrieve the data
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

    # multigen plot for tcells and tumors
    plot_settings = {}
    plot_agents_multigen(tcell_data, plot_settings, out_dir, tcell_id)
    plot_agents_multigen(tumor_data, plot_settings, out_dir, tumor_id)



# all of the experiments go here
simulation_experiments_library = {
    '1': simulation_1
}

# main experiment functions
def add_arguments():
    parser = argparse.ArgumentParser(description='tumor-tcell experiments')
    parser.add_argument(
        'experiment_id',
        type=str,
        choices=list(simulation_experiments_library.keys()),
        help='experiment id corresponds to ')
    return parser.parse_args()

def make_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def run_experiment():
    """ execute experiments """
    args = add_arguments()
    make_dir(EXPERIMENT_OUT_DIR)

    if args.experiment_id:
        # retrieve preset experiment
        experiment_id = str(args.experiment_id)
        experiment_type = simulation_experiments_library[experiment_id]
        control_out_dir = os.path.join(EXPERIMENT_OUT_DIR, experiment_id)
        make_dir(control_out_dir)

        if callable(experiment_type):
            experiment_type(control_out_dir)

    else:
        print('provide experiment number')


if __name__ == '__main__':
    run_experiment()
