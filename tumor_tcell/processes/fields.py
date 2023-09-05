"""
===============
Diffusion Field
===============

Diffuses and decays molecular concentrations in a 2D field.
"""

import os
import copy
import cv2
import numpy as np
from scipy import constants
# from scipy.ndimage import convolve

from vivarium.core.serialize import Quantity
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
CONCENTRATION_UNIT = units.ng / units.mL
LENGTH_UNIT = units.um
DIFFUSION_RATES = {
    'IFNg': 1.25e-3 * units.cm * units.cm / units.day,  # 1.25e-3 cm^2/day (Liao, 2014)
    'tumor_debris': 0.0864 * units.cm * units.cm / units.day,  # 10*10e-11 m^2/s (Krouglova, 2004)
}


class Fields(Process):
    """
    Diffusion and decay in 2-dimensional fields of molecules with agent exchange of molecules

    Arguments:
        parameters(dict): Accepts the following configuration keys:

        * **bounds** (list): size of the environment in micrometers, with ``[x, y]``.
        * **n_bins** (list): resolution of the environment by number of bins ``[i, j]``.
        * **depth** (Quantity): depth of the field.
        * **molecules** (list): the molecules, each will have its own field/
        * **default_diffusion_dt** (float): the time step with which to run diffusion.
            Must be less that the process time step.
        * **default_diffusion_rate** (float): the diffusion rate for all molecules that do not have a specific rate.
        * **diffusion** (dict): Specific diffusion rates for molecules with {'mol_id': rate}.
        * **decay** (dict): Specific decay rates for molecules with {'mol_id': rate}.
            If not provided, molecule does not decay.

    # TODO add recycling - 100-1000 molecules/cell/min #(Zhou, 2018)
    """

    name = NAME
    defaults = {
        # parameters for the lattice dimensions
        'bounds': [10 * units.um, 10 * units.um],
        'n_bins': [10, 10],
        'depth': 5000.0 * units.um,

        # molecules
        'molecules': ['IFNg'],

        # diffusion
        'default_diffusion_dt': 0.1,
        'default_diffusion_rate': 1e-1,
        # specific diffusion rates
        'diffusion': DIFFUSION_RATES,

        # specific decay rates
        'decay': {
            'IFNg': np.log(2)/(4.5*60*60),  # 7 hr half-life converted to exponential decay rate #(Kurzrock, 1985)
        },
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        # parameters
        self.molecule_ids = self.parameters['molecules']
        self.n_bins = self.parameters['n_bins']
        self.bounds = [
            b.to(LENGTH_UNIT).magnitude
            for b in self.parameters['bounds']]
        self.depth = self.parameters['depth']
        if isinstance(self.depth, Quantity):
            self.depth = self.depth.to(LENGTH_UNIT).magnitude


        # get diffusion rates
        diffusion_rate = self.parameters['default_diffusion_rate']
        bins_x = self.n_bins[0]
        bins_y = self.n_bins[1]
        length_x = self.bounds[0]
        length_y = self.bounds[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy

        # general diffusion rate
        self.diffusion_rate = diffusion_rate / dx2

        # diffusion rates for each individual molecules
        self.molecule_specific_diffusion = {
            mol_id: diff_rate.to(LENGTH_UNIT**2/units.s).magnitude/dx2
            for mol_id, diff_rate in self.parameters['diffusion'].items()}

        # get diffusion timestep
        diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * diffusion_rate * (dx ** 2 + dy ** 2))
        self.diffusion_dt = min(diffusion_dt, self.parameters['default_diffusion_dt'])

        # get bin volume, to convert between counts and concentration
        self.bin_volume = get_bin_volume(self.n_bins, self.bounds, self.depth)

    def initial_state(self, config=None):
        """get initial state of the fields

        Args:
            * config (dict): with optional keys "random" or "uniform".
                * "random" key maps to a maximum value for the field, which gets filled with values between [0, max].
                * "uniform" key maps to a value that will fill the entire field
        Returns:
            * fields (dict) with {mol_id: 2D np.array}
        """

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
        schema = {}

        # cells
        local_concentration_schema = {
            molecule: {
                '_default': 0.0}
            for molecule in self.parameters['molecules']}
        schema['cells'] = {
            '*': {
                'boundary': {
                    'location': {
                        '_default': [
                            0.5 * bound for bound in self.parameters['bounds']],
                     },
                    'external': local_concentration_schema
                }}}

        # fields
        fields_schema = {
            'fields': {
                field: {
                    '_default': self.ones_field(),
                    '_updater': 'nonnegative_accumulate',
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
            local_environment[mol_id] = field[bin_site]
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
        """ diffuse a single field """
        t = 0.0
        dt = min(timestep, self.diffusion_dt)
        diffusion_rate_dt = diffusion_rate * dt
        while t < timestep:
            result = cv2.filter2D(field, -1, LAPLACIAN_2D)
            # result = convolve(field, LAPLACIAN_2D, mode='reflect')
            field += diffusion_rate_dt * result
            t += dt
        return field

    def diffuse_fields(self, fields, timestep):
        """ diffuse fields in a fields dictionary """
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


def test_fields(config=None, initial=None, total_time=30):
    config = config or {'molecules': ['IFNg']}
    initial = initial or {'random': 5.0}
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

def main():
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    os.makedirs(out_dir, exist_ok=True)

    # run decay
    decay_config = {
        'molecules': ['IFNg'],
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
        'molecules': ['IFNg'],
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


if __name__ == '__main__':
    main()