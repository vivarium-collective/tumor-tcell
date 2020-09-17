'''
===============
Diffusion Field
===============
'''

import numpy as np
from scipy import constants

from vivarium.core.process import Process
from vivarium.library.units import units



NAME = 'diffusion_field'

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])
AVOGADRO = constants.N_A



class Fields(Process):
    '''
    Diffusion in 2-dimensional fields of molecules with agent exchange
    '''

    name = NAME
    defaults = {
        'time_step': 1,
        'molecules': ['glc'],
        'initial_state': {},
        'n_bins': [10, 10],
        'bounds': [10, 10],
        'depth': 3000.0,  # um
        'diffusion': 5e-1,
        'gradient': {},
    }

    def __init__(self, parameters=None):
        super(Fields, self).__init__(parameters)

    def ports_schema(self):
        local_concentration_schema = {
            molecule: {
                '_default': 0.0,
                '_updater': 'set'}
            for molecule in self.molecule_ids}

        schema = {'cells': {}}
        for agent_id, states in self.initial_agents.items():
            location = states['boundary'].get('location', [])
            schema['cells'][agent_id] = {
                'boundary': {
                    'location': {
                        '_value': location},
                }
            }
        glob_schema = {
            '*': {
                'boundary': {
                    'location': {
                        '_default': [0.5 * bound for bound in self.bounds],
                        '_updater': 'set'},
                    'external': local_concentration_schema}}}
        schema['cells'].update(glob_schema)

        # fields
        fields_schema = {
            'fields': {
                field: {
                    '_value': self.initial_state.get(field, self.ones_field()),
                    '_updater': 'accumulate',
                    '_emit': True,
                }
                for field in self.molecule_ids
            },
        }
        schema.update(fields_schema)

        # dimensions
        dimensions_schema = {
            'dimensions': {
                'bounds': {
                    '_value': self.parameters['bounds'],
                    '_updater': 'set',
                    '_emit': True,
                },
                'n_bins': {
                    '_value': self.parameters['n_bins'],
                    '_updater': 'set',
                    '_emit': True,
                },
                'depth': {
                    '_value': self.parameters['depth'],
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
        local_environments = self.get_local_environments(cells, fields)

        update = {'fields': delta_fields}
        if local_environments:
            update.update({'cells': local_environments})

        return update

    def count_to_concentration(self, count):
        return count_to_concentration(
            count * units.count, self.bin_volume * units.L
        ).to(units.mmol / units.L).magnitude

    def get_bin_site(self, location):
        return get_bin_site(location, self.n_bins, self.bounds)

    def get_single_local_environments(self, specs, fields):
        bin_site = self.get_bin_site(specs['location'])
        local_environment = {}
        for mol_id, field in fields.items():
            local_environment[mol_id] = field[bin_site]
        return local_environment

    def get_local_environments(self, cells, fields):
        local_environments = {}
        if cells:
            for agent_id, specs in cells.items():
                local_environments[agent_id] = {'boundary': {}}
                local_environments[agent_id]['boundary']['external'] = \
                    self.get_single_local_environments(specs['boundary'], fields)
        return local_environments

    def ones_field(self):
        return np.ones((self.n_bins[0], self.n_bins[1]), dtype=np.float64)

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
