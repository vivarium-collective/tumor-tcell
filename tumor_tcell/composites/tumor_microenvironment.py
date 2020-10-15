"""
======================
Tumor microenvironment
======================
"""

import os
import random

from vivarium.core.process import Generator
from vivarium.core.composition import (
    compartment_in_experiment,
    COMPOSITE_OUT_DIR,
)
from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units, remove_units

# processes
from tumor_tcell.processes.neighbors import Neighbors
from tumor_tcell.processes.fields import Fields
# from vivarium_cell.processes.diffusion_field import DiffusionField

# plots
from vivarium_cell.plots.multibody_physics import plot_snapshots

NAME = 'tumor_microenvironment'
DEFAULT_BOUNDS = [50 * units.um, 50 * units.um]


class TumorMicroEnvironment(Generator):
    """ Tumor micro-environment

    Models a spatial environment in which t-cells and tumors interact
    """
    name = NAME
    defaults = {
        'neighbors_multibody': {
            'bounds': DEFAULT_BOUNDS
        },
        'diffusion_field': {},
        '_schema': {},
    }

    def __init__(self, config=None):
        super(TumorMicroEnvironment, self).__init__(config)

    def initial_state(self, config=None):
        diffusion_field = Fields(self.config['diffusion_field'])
        return diffusion_field.initial_state(config)

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
        config['diffusion_field']['diffusion'] = diffusion
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
        config=None,
        n_agents=1,
        end_time=10
):
    if config is None:
        config = make_neighbors_config()
    # configure the compartment
    compartment = TumorMicroEnvironment(config)

    # set initial agent state
    initial_state = compartment.initial_state({'gradient': 'random'})
    if n_agents:
        agent_ids = [str(agent_id) for agent_id in range(n_agents)]
        body_config = {'agent_ids': agent_ids}
        if 'neighbors_multibody' in config and 'bounds' in config['neighbors_multibody']:
            body_config.update({'bounds': config['neighbors_multibody']['bounds']})
        initial_agents_state = agent_body_config(body_config)
        initial_state.update({'agents': initial_agents_state})

    # configure experiment
    experiment_settings = {
        'compartment': config,
        'initial_state': initial_state}
    experiment = compartment_in_experiment(
        compartment,
        experiment_settings)

    # run experiment
    experiment.update(end_time)
    data = experiment.emitter.get_data()

    # assert that the agent remains in the simulation until the end
    assert len(data[end_time]['agents']) == n_agents
    return data


def main():
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    bounds = [25 * units.um, 25 * units.um]
    config = make_neighbors_config(bounds=bounds)
    data = test_microenvironment(
        config=config,
        n_agents=1,
        end_time=40)

    # make snapshot plot
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_data = {
        'agents': agents,
        'fields': fields,
        'config': {'bounds': remove_units(bounds)}}
    plot_config = {
        'out_dir': out_dir,
        'agent_shape': 'circle',
        'filename': 'snapshots'}
    plot_snapshots(plot_data, plot_config)


if __name__ == '__main__':
    main()
