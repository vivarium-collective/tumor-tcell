"""
===========
Tumor Agent
===========
"""

import os

# core imports
from vivarium.core.process import Generator
from vivarium.core.composition import simulate_compartment_in_experiment
from vivarium.plots.agents_multigen import plot_agents_multigen

# processes
from vivarium.processes.meta_division import MetaDivision
from vivarium.processes.disintegrate import Disintegrate
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
        'death': {},
        '_schema': {},
    }

    def __init__(self, config):
        super(TumorAgent, self).__init__(config)

    def initial_state(self, config=None):
        process = TumorProcess()
        return process.initial_state(config)

    def generate_processes(self, config):
        daughter_path = config['daughter_path']
        agent_id = config['agent_id']
        meta_division_config = dict(
            {},
            daughter_path=daughter_path,
            agent_id=agent_id,
            compartment=self)

        return {
            'tumor': TumorProcess(config['tumor']),
            'division': MetaDivision(meta_division_config),
            # 'death': Disintegrate(),
        }

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']
        death_trigger_path = boundary_path + ('death',)
        return {
            'tumor': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'globals': boundary_path,
                'neighbors': ('neighbors',),
            },
            'division': {
                'global': boundary_path,
                'agents': agents_path,
            },
            # 'death': {
            #     'trigger': death_trigger_path,
            #     'agents': agents_path}
        }


# tests
def test_tumor_agent(total_time=1000):
    agent_id = '0'
    parameters = {'agent_id': agent_id}
    compartment = TumorAgent(parameters)

    # settings for simulation and plot
    settings = {
        'initial_state': compartment.initial_state(),
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
