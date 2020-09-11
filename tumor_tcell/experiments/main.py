"""
==============
Experiments
==============
"""

import os
import argparse

from vivarium.core.composition import agent_environment_experiment

from tumor_tcell.composites.tumor import TumorCompartment
from tumor_tcell.composites.t_cell import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment

# directories
from tumor_tcell import EXPERIMENT_OUT_DIR


def simulation_1():
    # configure the cells
    cell_config = [{
        'number': 1,
        'ids': ['tumor'],
        'type': TumorCompartment,
        'config': {},
        },
        {
        'number': 1,
        'ids': ['tcell'],
        'type': TCellAgent,
        'config': {},
        }]

    # configure the environment
    environment_config = {
        'type': TumorMicroEnvironment,
        'config': {},
    }

    # make the experiment using helper function agent_environment_experiment
    experiment = agent_environment_experiment(
        agents_config=cell_config,
        environment_config=environment_config
    )

    import ipdb; ipdb.set_trace()



# all of the experiments go here
simulation_experiments_library = {
    '1': simulation_1
}

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
    out_dir = os.path.join(EXPERIMENT_OUT_DIR)
    make_dir(out_dir)
    args = add_arguments()

    if args.experiment_id:
        # retrieve preset experiment
        experiment_id = str(args.experiment_id)
        experiment_type = simulation_experiments_library[experiment_id]
        control_out_dir = os.path.join(out_dir, experiment_id)
        make_dir(control_out_dir)

        if callable(experiment_type):
            data = experiment_type()

    else:
        print('provide experiment number')


if __name__ == '__main__':
    run_experiment()
