"""
============
T-cell Agent
============
"""

import os

# core imports
from vivarium.core.process import Generator
from vivarium.core.composition import simulate_compartment_in_experiment
from vivarium.plots.agents_multigen import plot_agents_multigen

# processes
from vivarium.processes.meta_division import MetaDivision
from vivarium.processes.disintegrate import Disintegrate
from tumor_tcell.processes.t_cell import TCellProcess
from tumor_tcell.processes.local_field import LocalField

# directories
from tumor_tcell import COMPOSITE_OUT_DIR


NAME = 'tcell_agent'

def daughter_ab(mother_id):
    return [
        str(mother_id) + "A",
        str(mother_id) + "B"
    ]

class TCellAgent(Generator):

    name = NAME
    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        'field_path': ('..', '..', 'fields',),
        'dimensions_path': ('..', '..', 'dimensions',),
        'tcell': {},
        'death': {},
        '_schema': {
            'death': {
                'trigger': {
                    '_emit': False
                }
            },
            't_cell': {
                'globals': {
                    'death': {
                        '_emit': False,
                    },
                    'divide': {
                        '_emit': False,
                    },
                    'PD1n_divide_count': {
                        '_emit': False,
                    },
                    'PD1p_divide_count': {
                        '_emit': False,
                    },
                },
                'internal': {
                    'cell_state_count': {
                        '_emit': False,
                    }
                },
            },
        },
    }

    def __init__(self, config):
        super(TCellAgent, self).__init__(config)

    def initial_state(self, config=None):
        process = TCellProcess()
        return process.initial_state(config)

    def generate_processes(self, config):
        daughter_path = config['daughter_path']
        agent_id = config['agent_id']

        # division config
        meta_division_config = dict(
            {},
            daughter_ids_function = daughter_ab,
            daughter_path=daughter_path,
            agent_id=agent_id,
            compartment=self)

        # death config
        death_config = {
            'agent_id': agent_id}

        return {
            't_cell': TCellProcess(config['tcell']),
            'division': MetaDivision(meta_division_config),
            'death': Disintegrate(death_config),
            'local_field': LocalField(),
        }

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']
        field_path = config['field_path']
        dimensions_path = config['dimensions_path']
        death_trigger_path = boundary_path + ('death',)

        return {
            't_cell': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'globals': boundary_path,
                'neighbors': ('neighbors',),
            },
            'local_field': {
                'exchanges': boundary_path + ('external',),
                'location': boundary_path + ('location', ),
                'fields': field_path,
                'dimensions': dimensions_path,
            },
            'division': {
                'global': boundary_path,
                'agents': agents_path,
            },
            'death': {
                'trigger': death_trigger_path,
                'agents': agents_path}
        }


# tests
def test_tcell_agent(total_time=1000):
    agent_id = '0'
    parameters = {'agent_id': agent_id}
    compartment = TCellAgent(parameters)

    # settings for simulation and plot
    settings = {
        'initial_state': compartment.initial_state(),
        'outer_path': ('agents', agent_id),
        'return_raw_data': True,
        'timestep': 10,
        'total_time': total_time}
    return simulate_compartment_in_experiment(compartment, settings)

def run_compartment(out_dir='out'):
    data = test_tcell_agent(total_time=100000)
    plot_settings = {}
    plot_agents_multigen(data, plot_settings, out_dir)


if __name__ == '__main__':
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    run_compartment(out_dir=out_dir)
