"""
====================
Local Fields Deriver
====================

This deriver helps translate a large 2D environment into an agent's local environment,
so that it can access the local concentrations without any knowledge of the larger
environmental state.
"""

import numpy as np
from scipy import constants

from vivarium.core.process import Deriver
from vivarium.library.units import units, remove_units

from vivarium_multibody.library.lattice_utils import (
    get_bin_site,
    get_bin_volume,
    count_to_concentration,
)

CONCENTRATION_UNIT = units.ng / units.mL  # alternative (units.mmol / units.L) concentration would not use molecular_weight
LENGTH_UNIT = units.um
UNITLESS_AVOGADRO = constants.N_A


class LocalField(Deriver):

    name = 'local_field'
    defaults = {
        'molecular_weight': {
            'IFNg': 17000 * units.g / units.mol,
            'tumor_debris': 29000 * units.g / units.mol, #Assumed HGMB1 (Apetoh, 2007)
        },
        'concentration_unit': CONCENTRATION_UNIT
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        self.conc_conversion = {}
        for mol_id, mw in self.parameters['molecular_weight'].items():
            self.conc_conversion[mol_id] = (units.mol / units.L * mw).to(
                self.parameters['concentration_unit']).magnitude

    def ports_schema(self):
         return {
            'exchanges': {
                '*': {
                    '_default': 0,  # counts!
                }
            },
            'location': {
                '_default': [0.5 * LENGTH_UNIT, 0.5 * LENGTH_UNIT]
            },
            'fields': {
                '*': {
                    '_default': np.ones(1),
                    #'_updater': 'accumulate',
                }
            },
            'dimensions': {
                'bounds': {
                    '_default': [1, 1],
                },
                'n_bins': {
                    '_default': [1, 1],
                },
                'depth': {
                    '_default': 1,
                },
            }
        }


    def next_update(self, timestep, states):
        if not states:
            return {}
        location = remove_units(states['location'])
        n_bins = states['dimensions']['n_bins']
        bounds = states['dimensions']['bounds']
        depth = states['dimensions']['depth']
        exchanges = states['exchanges']

        # get bin
        bin_site = get_bin_site(location, n_bins, bounds)
        bin_volume_liters = get_bin_volume(n_bins, bounds, depth)

        # apply exchanges
        delta_fields = {}
        reset_exchanges = {}
        for mol_id, counts in exchanges.items():
            delta_fields[mol_id] = np.zeros(
                (n_bins[0], n_bins[1]), dtype=np.float64)
            concentration = counts / (bin_volume_liters * UNITLESS_AVOGADRO) * self.conc_conversion[mol_id]
            delta_fields[mol_id][bin_site[0], bin_site[1]] += concentration
            reset_exchanges[mol_id] = {
                '_value': 0,
                '_updater': 'set'}

        return {
            'exchanges': reset_exchanges,
            'fields': delta_fields}


def test_local_fields():
    mol_name = 'A'
    parameters = {'molecular_weight': {mol_name: 1e3 * units.g / units.mol}}
    local_fields_process = LocalField(parameters)

    bounds = [5, 5]
    n_bins = [3, 3]
    initial_state = {
        'exchanges': {
            mol_name: 20
        },
        'location': [0.5, 0.5],
        'fields': {
            mol_name: np.ones((n_bins[0], n_bins[1]), dtype=np.float64)
        },
        'dimensions': {
            'bounds': bounds,
            'n_bins': n_bins,
            'depth': 1,
        }
    }

    output = local_fields_process.next_update(0, initial_state)


if __name__ == '__main__':
    test_local_fields()
