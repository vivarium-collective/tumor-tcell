'''
===============
Diffusion Field
===============
'''

import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from vivarium.core.process import Process
from vivarium.library.units import units
from vivarium_cell.library.lattice_utils import (
    count_to_concentration,
    get_bin_site,
    get_bin_volume,
)


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
        'molecules': ['molecule'],
        'initial_state': {},
        'n_bins': [10, 10],
        'bounds': [10 * units.um, 10 * units.um],
        'depth': 3000.0,  # um
        'diffusion': 5e-1,
        'gradient': {},

        #TODO - Add the diffusion and degradation of IFNg
        'IFNg_diffusion_coeff': 1.25 * 10**-3, #cm^2/day #(Liao, 2014)
        'IFNg_degradation_coeff': 2.16,  # 1/day #(Liao, 2014)
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

        # diffusion
        diffusion = self.parameters['diffusion']
        bins_x = self.n_bins[0]
        bins_y = self.n_bins[1]
        length_x = self.bounds[0]
        length_y = self.bounds[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        self.diffusion = diffusion / dx2
        self.diffusion_dt = 0.01
        # self.diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * self.diffusion * (dx ** 2 + dy ** 2))

        # volume, to convert between counts and concentration
        self.bin_volume = get_bin_volume(self.n_bins, self.bounds, self.depth)

    def initial_state(self, config=None):
        if config is None:
            config = {}
        gradient = config.get('gradient', 'ones')
        if gradient == 'random':
            return {
                'fields': {
                    field: self.random_field()
                    for field in self.parameters['molecules']}}
        else:
            return {
                'fields': {
                    field: self.ones_field()
                    for field in self.parameters['molecules']}}

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

        # get neighbors
        # secrete
        # t-cells secrete a count of IFNg, needs to convert to concentration for tumors
        # TODO -- IFNg gets secreted from t-cells in the environment based on tumor neighbors (in fields process)

        # get each agent's local environment
        local_environments = self.set_local_environments(cells, fields)

        update = {'fields': delta_fields}
        if local_environments:
            update.update({'cells': local_environments})

        return update

    def count_to_concentration(self, count):
        return count_to_concentration(
            count * units.count, self.bin_volume * units.L
        ).to(units.mmol / units.L).magnitude

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