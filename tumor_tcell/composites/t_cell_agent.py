"""
============
T-cell Agent
============

This composite model combines the t-cell process with local_field, division, and death.

This composite can be run on its own from the command line.

    $ python tumor_tcell/composites/t_cell_agent.py

"""

import os

# core imports
from vivarium.core.composer import Composer, Composite
from vivarium.core.engine import Engine
from vivarium.plots.agents_multigen import plot_agents_multigen

# processes
from vivarium.processes.meta_division import MetaDivision
from vivarium.processes.timeline import TimelineProcess
from vivarium.processes.remove import Remove
from tumor_tcell.processes.t_cell import TCellProcess, TIMESTEP
from tumor_tcell.processes.local_field import LocalField

# directories/libraries
from tumor_tcell.library.phylogeny import daughter_ab
from tumor_tcell import COMPOSITE_OUT_DIR


NAME = 'tcell_agent'


class TCellAgent(Composer):

    name = NAME
    defaults = {
        'reuse_processes': False,
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
                        '_emit': True,
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
        super().__init__(config)
        self.processes_initialized = False

    def initialize_processes(self, config):
        self.tcell_process = TCellProcess(config['tcell'])
        self.local_field = LocalField()

        if self.config['reuse_processes']:
            self.processes_initialized = True

    def generate_processes(self, config):
        if not self.processes_initialized:
            self.initialize_processes(config)

        daughter_path = config['daughter_path']
        agent_id = config['agent_id']

        # division config
        meta_division_config = dict(
            {},
            daughter_ids_function=daughter_ab,
            daughter_path=daughter_path,
            agent_id=agent_id,
            composer=self)

        # death config
        death_config = {
            'agent_id': agent_id}

        return {
            't_cell': self.tcell_process,
            'local_field': self.local_field,
            'division': MetaDivision(meta_division_config),
            'death': Remove(death_config),
        }

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']
        field_path = config['field_path']
        dimensions_path = config['dimensions_path']

        return {
            't_cell': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'globals': boundary_path,
                'neighbors': ('neighbors',),
            },
            'local_field': {
                'exchanges': boundary_path + ('exchange',),
                'location': boundary_path + ('location',),
                'fields': field_path,
                'dimensions': dimensions_path,
            },
            'division': {
                'global': boundary_path,
                'agents': agents_path,
            },
            'death': {
                'trigger': boundary_path + ('death',),
                'agents': agents_path,
            },
        }


# tests
def test_tcell_agent(
        total_time=1000,
        agent_ids=['0'],
        agent_timeline=None,
        initial_agent_state='PD1n',
):
    composite = Composite()
    for agent_id in agent_ids:
        parameters = {
            'agent_id': agent_id,
            '_schema': {
                't_cell': {
                    'internal': {
                        'cell_state': {'_emit': False},
                        'velocity_timer': {'_emit': False}},
                    'neighbors': {
                        'present': {
                            'TCR': {'_emit': False}},
                        'transfer': {
                            'cytotoxic_packets': {'_emit': False}}},
                    'boundary': {
                        'velocity': {'_emit': False},
                        'external': {
                            'IFNg': {'_emit': False}}},
                    'globals': {
                        'PD1n_divide_count': {'_emit': True},
                        'PD1p_divide_count': {'_emit': True}}}}}

        composer = TCellAgent(parameters)
        agent = composer.generate(path=('agents', agent_id))
        composite.merge(composite=agent)

    if agent_timeline:
        timeline = []
        for agent_id in agent_ids:
            individual_timeline = []
            for (t, perturb) in agent_timeline:
                agent_perturb = {}
                for path, magnitude in perturb.items():
                    agent_path = ('agents', agent_id) + path
                    agent_perturb[agent_path] = magnitude
                agent_event = (t, agent_perturb)
                individual_timeline.append(agent_event)
            timeline.extend(individual_timeline)
            
        timeline_process = TimelineProcess({'timeline': timeline, 'time_step': TIMESTEP})
        composite.merge(composite=timeline_process.generate())

    # settings for simulation and plot
    initial = composite.initial_state()
    for agent_id in agent_ids:
        initial['agents'][agent_id]['internal']['cell_state'] = initial_agent_state  # set an initial state

    # make the experiment
    experiment = Engine(**{
        'processes': composite['processes'],
        'topology': composite['topology'],
        'initial_state': initial})
    experiment.update(total_time)

    output = experiment.emitter.get_data()

    # convert time to hours, and add agents key back in if all agents have died
    times = list(output.keys())
    for t in times:
        if 'agents' not in output[t]:
            output[t]['agents'] = {}
        hr = t/3600
        output[hr] = output.pop(t)

    return output

def run_agent(out_dir='out'):
    agent_ids = ['0', '1']
    agent_timeline = [
        (500, {
            ('neighbors', 'accept', 'PDL1'): 5e5,
            ('neighbors', 'accept', 'MHCI'): 5e5,
        }),
    ]
    data = test_tcell_agent(
        total_time=100000,
        agent_ids=agent_ids,
        agent_timeline=agent_timeline)
    plot_settings = {
        'time_display': '(hr)'
    }
    plot_agents_multigen(data, plot_settings, out_dir)


if __name__ == '__main__':
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    run_agent(out_dir=out_dir)
