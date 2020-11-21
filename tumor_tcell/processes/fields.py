'''
===============
Diffusion Field
===============
'''

import os

import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from vivarium.core.process import Process
from vivarium.core.composition import simulate_process
from vivarium.library.units import units, remove_units
from vivarium_cell.library.lattice_utils import (
    count_to_concentration,
    get_bin_site,
    get_bin_volume,
)

# plotting
from vivarium_cell.plots.multibody_physics import plot_snapshots

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'fields'

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])
AVOGADRO = constants.N_A

LENGTH_UNIT = units.um
CONCENTRATION_UNIT = 1  # TODO (ERAN) set value -- units.ng / units.mL


class Fields(Process):
    '''
    Diffusion in 2-dimensional fields of molecules with agent exchange
    '''

    name = NAME
    defaults = {
        'time_step': 1,
        'molecules': ['IFNg'],
        'initial_state': {},
        'n_bins': [10, 10],
        'bounds': [10 * units.um, 10 * units.um],
        'depth': 3000.0,  # um
        'default_diffusion_dt': 0.01,
        'default_diffusion_rate': 5e-1,
        'gradient': {},
        'diffusion': {
            'IFNg': 1.25e-3,  # cm^2/day #(Liao, 2014)
        },
        'degradation': {
            'IFNg': 2.16,  # 1/day #(Liao, 2014)
        },
        #If we want to add recycling - 100-1000 molecules/cell/min #(Zhou, 2018)
    }

    def __init__(self, parameters=None):
        super(Fields, self).__init__(parameters)

        # initial state
        self.molecule_ids = self.parameters['molecules']

        # parameters
        self.n_bins = self.parameters['n_bins']
        self.bounds = [b.to(LENGTH_UNIT).magnitude for b in self.parameters['bounds']]
        self.depth = self.parameters['depth']

        # get diffusion rates
        # TODO (ERAN) -- make molecule-specific diffusion rate dictionary
        diffusion_rate = self.parameters['default_diffusion_rate']
        bins_x = self.n_bins[0]
        bins_y = self.n_bins[1]
        length_x = self.bounds[0]
        length_y = self.bounds[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        self.diffusion = diffusion_rate / dx2

        # get diffusion timestep
        diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * diffusion_rate * (dx ** 2 + dy ** 2))
        self.diffusion_dt = min(diffusion_dt, self.parameters['default_diffusion_dt'])

        # volume, to convert between counts and concentration
        self.bin_volume = get_bin_volume(self.n_bins, self.bounds, self.depth)

    def initial_state(self, config=None):
        if config is None:
            config = {}
        if 'random' in config:
            max = config.get('random', 1)
            fields = {
                field: max * self.random_field()
                for field in self.parameters['molecules']}
        elif 'uniform' in config:
            fields = {
                field: config['uniform'] * self.ones_field()
                for field in self.parameters['molecules']}
        else:
            fields = {
                field: self.ones_field()
                for field in self.parameters['molecules']}
        return {
            'fields': fields,
            'cells': {},
        }

    def ports_schema(self):
        local_concentration_schema = {
            molecule: {
                '_default': 0.0 * CONCENTRATION_UNIT}
            for molecule in self.parameters['molecules']}

        schema = {}
        schema['cells'] = {
            '*': {
                'boundary': {
                    'location': {
                        '_default': [0.5 * bound for bound in self.bounds],
                        '_updater': 'set'},
                    'external': local_concentration_schema
                }}}

        # fields
        fields_schema = {
            'fields': {
                field: {
                    # '_value': self.initial_state.get(field, self.ones_field()),
                    '_default': self.ones_field(),
                    '_updater': 'accumulate',
                    '_emit': True,
                }
                for field in self.parameters['molecules']
            },
        }
        schema.update(fields_schema)

        # dimensions
        dimensions_schema = {
            'dimensions': {
                'bounds': {
                    '_value': self.bounds,
                    '_updater': 'set',
                    '_emit': True,
                },
                'n_bins': {
                    '_value': self.n_bins,
                    '_updater': 'set',
                    '_emit': True,
                },
                'depth': {
                    '_value': self.depth,
                    '_updater': 'set',
                    '_emit': True,
                }
            },
        }
        schema.update(dimensions_schema)
        return schema

    def next_update(self, timestep, states):
        fields = states['fields']
        cells = states['cells']

        # diffuse field
        delta_fields = self.diffuse(fields, timestep)

        # TODO -- add degradation

        # get each agent's local environment
        local_environments = self.set_local_environments(cells, fields)

        update = {'fields': delta_fields}
        if local_environments:
            update.update({'cells': local_environments})

        return update

    def get_bin_site(self, location):
        return get_bin_site(
            [l.to(LENGTH_UNIT).magnitude for l in location],
            self.n_bins,
            self.bounds)

    def get_single_local_environments(self, specs, fields):
        bin_site = self.get_bin_site(specs['location'])
        local_environment = {}
        for mol_id, field in fields.items():
            local_environment[mol_id] = field[bin_site] * CONCENTRATION_UNIT
        return local_environment

    def set_local_environments(self, cells, fields):
        local_environments = {}
        if cells:
            for agent_id, specs in cells.items():
                local_environments[agent_id] = {'boundary': {'external': {}}}
                cell_environment = self.get_single_local_environments(specs['boundary'], fields)
                local_environments[agent_id]['boundary']['external'] = {
                    mol_id: {
                        '_value': value,
                        '_updater': 'set'  # this overrides the default updater
                    } for mol_id, value in cell_environment.items()
                }
        return local_environments

    def ones_field(self):
        return np.ones((
            self.parameters['n_bins'][0],
            self.parameters['n_bins'][1]),
            dtype=np.float64)

    def random_field(self):
        return np.random.rand(
            self.parameters['n_bins'][0],
            self.parameters['n_bins'][1])

    # diffusion functions
    def diffusion_delta(self, field, timestep):
        ''' calculate concentration changes cause by diffusion'''
        field_new = field.copy()
        t = 0.0
        dt = min(timestep, self.diffusion_dt)
        while t < timestep:
            field_new += self.diffusion * dt * convolve(field_new, LAPLACIAN_2D, mode='reflect')
            t += dt

        return field_new - field

    def diffuse(self, fields, timestep):
        delta_fields = {}
        for mol_id, field in fields.items():
            # run diffusion if molecule field is not uniform
            if len(set(field.flatten())) != 1:
                delta = self.diffusion_delta(field, timestep)
            else:
                delta = np.zeros_like(field)
            delta_fields[mol_id] = delta

        return delta_fields

    def degrade(self, fields, timestep):
        import ipdb; ipdb.set_trace()
        return {}

def test_fields(config={}):

    # initialize process
    fields = Fields(config)

    # get initial state from process
    initial_state = fields.initial_state({'random': 5.0})

    # run experiment
    settings = {
        'return_raw_data': True,
        'initial_state': initial_state,
        'total_time': 10,
        'timestep': 1}
    return simulate_process(fields, settings)

def plot_fields(data, config, out_dir='out', filename='fields'):
    fields = {
        time: time_data['fields']
        for time, time_data in data.items()}
    snapshots_data = {
        'fields': fields,
        'config': config}
    plot_config = {
        'out_dir': out_dir,
        'filename': filename}
    plot_snapshots(snapshots_data, plot_config)


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    os.makedirs(out_dir, exist_ok=True)

    config = {
        'n_bins': [10, 10],
        'bounds': [10 * units.um, 10 * units.um]}
    data = test_fields(config)

    # plot
    plot_fields(
        data,
        remove_units(config),
        out_dir, 'test_field')
