"""
======================
Tumor microenvironment
======================

Composite model that combines neighbors and fields processes to make a lattice environment.

This composite can be run on its own from the command line.

    $ python tumor_tcell/composites/tumor_microenvironment.py

"""

import os
import random
import numpy as np

from vivarium.core.composer import Composer
from vivarium.core.engine import Engine
from vivarium.core.composition import COMPOSITE_OUT_DIR
from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units, remove_units

# processes
from tumor_tcell.processes.neighbors import Neighbors
from tumor_tcell.processes.fields import Fields
from tumor_tcell.processes.lymph_node import LymphNode

# plots
from tumor_tcell.plots.snapshots import plot_snapshots, format_snapshot_data

NAME = 'tumor_microenvironment'
DEFAULT_BOUNDS = [50 * units.um, 50 * units.um]


class TumorMicroEnvironment(Composer):
    """ Tumor micro-environment

    Models a spatial environment in which t-cells and tumors interact
    """
    name = NAME
    defaults = {
        'neighbors_multibody': {
            'name': 'neighbors_multibody',
            'bounds': DEFAULT_BOUNDS,
            'jitter_force': 0,
        },
        'diffusion_field': {
            'name': 'diffusion_field',
        },
        '_schema': {},
    }

    def __init__(self, config=None):
        super().__init__(config)

    def generate_processes(self, config):

        # initialize processes
        neighbors_multibody = Neighbors(config['neighbors_multibody'])
        diffusion_field = Fields(config['diffusion_field'])

        # make dictionary of processes
        return {
            'neighbors_multibody': neighbors_multibody,
            'diffusion_field': diffusion_field,
        }

    def generate_topology(self, config):
        return {
            'neighbors_multibody': {
                'cells': ('agents',)
            },
            'diffusion_field': {
                'cells': ('agents',),
                'fields': ('fields',),
                'dimensions': ('dimensions',),
            },
        }


class TumorAndLymphNodeEnvironment(Composer):
    """This Composer will have a 2D environment like TumorMicroEnvironment,
    and also a 2nd nonspatial lymph node environment that cells can migrate to"""
    defaults = {
        'tumor_env_id': 'tumor_environment',
        'ln_id': 'lymph_node',
        'neighbors_multibody': {
            'name': 'neighbors_multibody',
            'bounds': DEFAULT_BOUNDS,
            'jitter_force': 0,
        },
        'diffusion_field': {
            'name': 'diffusion_field',
        },
        'lymph_node': {
            'tumor_env_bounds': DEFAULT_BOUNDS,
        },
        '_schema': {},
    }

    def __init__(self, config=None):
        super().__init__(config)

    def generate_processes(self, config):
        # initialize processes
        neighbors_multibody = Neighbors(config['neighbors_multibody'])
        diffusion_field = Fields(config['diffusion_field'])
        lymph_node_transfer_process = LymphNode(config['lymph_node'])

        # make dictionary of processes
        return {
            config['tumor_env_id']: {
                'neighbors_multibody': neighbors_multibody,
                'diffusion_field': diffusion_field,
            },
            'lymph_node_transfer': lymph_node_transfer_process
        }

    def generate_topology(self, config):
        return {
            config['tumor_env_id']: {
                'neighbors_multibody': {
                    'cells': ('agents',)
                },
                'diffusion_field': {
                    'cells': ('agents',),
                    'fields': ('fields',),
                    'dimensions': ('dimensions',),
                },
            },
            'lymph_node_transfer': {
                'cells': (config['tumor_env_id'],),
                'lymph_node': (config['ln_id'],),
                'in_transit': ('in_transit',),
            }
        }


def make_neighbors_config(
        time_step=None,
        jitter_force=None,
        bounds=None,
        n_bins=None,
        depth=None,
        concentrations=None,
        molecules=None,
        diffusion=None,
        keep_fields_emit=None,
        set_config=None,
        parallel=None,
):
    config = {'neighbors_multibody': {}, 'diffusion_field': {}}

    if time_step:
        config['neighbors_multibody']['time_step'] = time_step
        config['diffusion_field']['time_step'] = time_step
    if bounds:
        config['neighbors_multibody']['bounds'] = bounds
        config['diffusion_field']['bounds'] = bounds
        config['diffusion_field']['n_bins'] = remove_units(bounds)
    if n_bins:
        config['diffusion_field']['n_bins'] = n_bins
    if jitter_force:
        config['neighbors_multibody']['jitter_force'] = jitter_force
    if depth:
        config['diffusion_field']['depth'] = depth
    if diffusion:
        config['diffusion_field']['default_diffusion_rate'] = diffusion
    if concentrations:
        config['diffusion_field']['gradient'] = {
            'type': 'uniform',
            'molecules': concentrations}
        molecules = list(concentrations.keys())
        config['diffusion_field']['molecules'] = molecules
    elif molecules:
        # molecules are a list, assume uniform concentrations of 1
        config['diffusion_field']['molecules'] = molecules
    if keep_fields_emit:
        # by default no fields are emitted
        config['diffusion_field']['_schema'] = {
            'fields': {
                field_id: {
                    '_emit': False}
                for field_id in molecules
                if field_id not in keep_fields_emit}}
    if parallel:
        config['diffusion_field']['_parallel'] = True
        config['neighbors_multibody']['_parallel'] = True
    if set_config:
        config = deep_merge(config, set_config)

    return config

def single_agent_config(config):
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    location = config.get('location')
    if location:
        location = [loc * bounds[n] for n, loc in enumerate(location)]
    else:
        location = [random.uniform(0, b) for b in bounds]

    return {'boundary': {
        'location': location,
        'diameter': 1 * units.um,
        'mass': 1339 * units.fg,
    }}

def agent_body_config(config):
    agent_ids = config['agent_ids']
    agent_config = {
        agent_id: single_agent_config(config)
        for agent_id in agent_ids}
    return agent_config

def test_microenvironment(
        n_agents=1,
        bounds=[25 * units.um, 25 * units.um],
        n_bins=[25, 25],
        jitter_force=1e-6,
        diffusion=1e-1,
        end_time=10
):

    # configure the compartment
    config = make_neighbors_config(
        bounds=bounds,
        n_bins=n_bins,
        jitter_force=jitter_force,
        diffusion=diffusion)
    compartment = TumorMicroEnvironment(config)

    # set initial agent state
    # initial_state = compartment.initial_state({'diffusion_field': {'random': 1}})
    initial_field = np.zeros((n_bins[0], n_bins[1]))
    initial_field[:, -1] = 100
    initial_state = {
        'fields': {'IFNg': initial_field},
        'agents': {}}
    if n_agents:
        agent_ids = [str(agent_id) for agent_id in range(n_agents)]
        body_config = {'agent_ids': agent_ids}
        if 'neighbors_multibody' in config and 'bounds' in config['neighbors_multibody']:
            body_config.update({'bounds': config['neighbors_multibody']['bounds']})
        initial_agents_state = agent_body_config(body_config)
        initial_state.update({'agents': initial_agents_state})

    # configure experiment
    composite = compartment.generate()
    experiment = Engine(
        composite=composite,
        initial_state=initial_state)

    # run experiment
    experiment.update(end_time)
    data = experiment.emitter.get_data_deserialized()

    # assert that the agent remains in the simulation until the end
    assert len(data[end_time]['agents']) == n_agents
    return data


def main():
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # configure
    bounds = [10 * units.um, 10 * units.um]
    n_bins = [10, 10]

    # run the simulation
    data = test_microenvironment(
        bounds=bounds,
        n_bins=n_bins,
        jitter_force=1e-4,
        diffusion=1e-9,
        n_agents=1,
        end_time=6)

    # snapshot plot
    agents, fields = format_snapshot_data(data)
    plot_snapshots(
        agents=remove_units(agents),
        fields=fields,
        bounds=remove_units(bounds),
        n_snapshots=5,
        out_dir=out_dir,
        filename='snapshots'
    )


if __name__ == '__main__':
    main()
