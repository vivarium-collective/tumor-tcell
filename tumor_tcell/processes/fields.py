'''
===============
Diffusion Field
===============
'''

import os
import copy

import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from vivarium.core.process import Process
from vivarium.core.composition import simulate_process
from vivarium.library.units import units, remove_units
from vivarium_multibody.library.lattice_utils import (
    get_bin_site,
    get_bin_volume,
)

# plotting
from tumor_tcell.plots.snapshots import plot_snapshots

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
        'depth': 5000.0,  # um
        'default_diffusion_dt': 0.1,
        'default_diffusion_rate': 1e-1,
        'gradient': {},
        'diffusion': {
            'IFNg': 1.25e-3 * units.cm * units.cm / units.day,  # 1.25e-3 cm^2/day #(Liao, 2014)
        },
        'decay': {
            'IFNg': np.log(2)/(4.5*60*60),  # 4.5 hr half-life converted to exponential decay rate #(Kurzrock, 1985)
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
        self.diffusion_rate = diffusion_rate / dx2
        self.molecule_specific_diffusion = {
            mol_id: diff_rate.to(LENGTH_UNIT**2/units.s).magnitude/dx2
            for mol_id, diff_rate in self.parameters['diffusion'].items()}

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
                        # '_default': [0.5 * bound for bound in self.bounds],
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

        # degrade and diffuse
        # TODO -- use set update, so you don't need to do the copy and get the delta
        fields_new = copy.deepcopy(fields)
        fields_new = self.degrade_fields(fields_new, timestep)
        fields_new = self.diffuse_fields(fields_new, timestep)

        # get delta_fields
        delta_fields = {
            mol_id: fields_new[mol_id] - field
            for mol_id, field in fields.items()}

        # get each agent's local environment
        local_environments = self.set_local_environments(cells, fields_new)

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

    def diffuse(self, field, timestep, diffusion_rate):
        ''' diffuse a single field '''
        t = 0.0
        dt = min(timestep, self.diffusion_dt)
        while t < timestep:
            field += diffusion_rate * dt * convolve(field, LAPLACIAN_2D, mode='reflect')
            t += dt
        return field

    def diffuse_fields(self, fields, timestep):
        ''' diffuse fields in a fields dictionary '''
        for mol_id, field in fields.items():
            diffusion_rate = self.molecule_specific_diffusion.get(mol_id, self.diffusion_rate)
            # run diffusion if molecule field is not uniform
            if len(set(field.flatten())) != 1:
                fields[mol_id] = self.diffuse(field, timestep, diffusion_rate)
        return fields

    def degrade_fields(self, fields, timestep):
        """
        Note: this only applies if the molecule has a rate in parameters['decay']
        """
        for mol_id, field in fields.items():
            if mol_id in self.parameters['decay']:
                decay_rate = self.parameters['decay'][mol_id]
                degraded = field * np.exp(-decay_rate * timestep)
                fields[mol_id] = degraded
        return fields

def test_fields(config={}, initial={'random': 5.0}, total_time=30):
    # initialize process
    fields = Fields(config)

    # get initial state from process
    initial_state = fields.initial_state(initial)

    # run experiment
    settings = {
        'return_raw_data': True,
        'initial_state': initial_state,
        'total_time': total_time,
    }
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

    # run decay
    decay_config = {
        'time_step': 60,
        'n_bins': [10, 10],
        'bounds': [10 * units.um, 10 * units.um]}
    decay_data = test_fields(
        config=decay_config,
        initial={'uniform': 1.0},
        total_time=4.5*60*60)
    plot_fields(
        decay_data,
        remove_units(decay_config),
        out_dir, 'test_decay')

    # run diffuse
    diffuse_config = {
        'time_step': 60,
        'n_bins': [10, 10],
        'bounds': [10 * units.um, 10 * units.um]}
    diffuse_data = test_fields(
        config=diffuse_config,
        initial={'random': 1.0},
        total_time=4.5*60*60)
    plot_fields(
        diffuse_data,
        remove_units(diffuse_config),
        out_dir, 'test_diffuse')