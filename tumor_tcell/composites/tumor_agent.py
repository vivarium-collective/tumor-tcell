"""
===========
Tumor Agent
===========
"""

import os

from vivarium.core.composition import simulate_compartment_in_experiment
from vivarium.plots.agents_multigen import plot_agents_multigen

# processes
from vivarium.core.process import Generator
from vivarium.processes.meta_division import MetaDivision
from tumor_tcell.processes.tumor import TumorProcess

# directories
from tumor_tcell import COMPOSITE_OUT_DIR

NAME = 'tumor_agent'


class TumorAgent(Generator):

    name = NAME
    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        'tumor': {},
        'divide': True,
        '_schema': {},
    }

    def __init__(self, config):
        super(TumorAgent, self).__init__(config)

    def generate_processes(self, config):

        # initialize processes
        tumor = TumorProcess(config['tumor'])

        # make dictionary of processes
        processes = {
            'tumor': tumor,
        }

        # if divide set to true, add meta-division processes
        if config['divide']:
            daughter_path = config['daughter_path']
            agent_id = config['agent_id']
            meta_division_config = dict(
                {},
                daughter_path=daughter_path,
                agent_id=agent_id,
                compartment=self)
            meta_division = MetaDivision(meta_division_config)
            processes['division'] = meta_division

        return processes


    def generate_topology(self, config):

        # retrieve paths
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']

        # make topology by mapping ports
        topology = {
            'tumor': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'globals': boundary_path,
                'neighbors': ('neighbors',),
            },
        }
        if config['divide']:
            topology.update({
                'division': {
                    'global': boundary_path,
                    'agents': agents_path,
                }})
        return topology


# tests
def test_tumor_agent(total_time=1000):
    agent_id = '0'
    parameters = {'agent_id': agent_id}
    compartment = TumorAgent(parameters)

    # settings for simulation and plot
    settings = {
        'outer_path': ('agents', agent_id),
        'return_raw_data': True,
        'timestep': 10,
        'total_time': total_time}
    return simulate_compartment_in_experiment(compartment, settings)

def run_compartment(out_dir='out'):
    data = test_tumor_agent(total_time=4000)
    plot_settings = {}
    plot_agents_multigen(data, plot_settings, out_dir)


if __name__ == '__main__':
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    run_compartment(out_dir=out_dir)
