"""
==============
Dendritic Cell Process
==============
"""

import random
from scipy import constants

from vivarium.core.process import Process
from vivarium.core.composer import Composite
from vivarium.core.engine import Engine
from vivarium.plots.simulation_output import plot_simulation_output
from vivarium.library.units import units
from vivarium.processes.timeline import TimelineProcess

from tumor_tcell.processes.fields import DIFFUSION_RATES  # , CONCENTRATION_UNIT
from tumor_tcell.processes.tumor import get_probability_timestep

TIMESTEP = 60  # seconds
CONCENTRATION_UNIT = units.ng / units.mL  # TODO: units.ng / units.mL
AVOGADRO = constants.N_A
PI = constants.pi


class DendriticCellProcess(Process):
    """DendriticCellProcess

    References:
        *
    """
    defaults = {
        'mass': 2.0 * units.ng,  # TODO
        'diameter': 10.0 * units.um,  # * units.um, (Morefield, 2005)
        'velocity': 3.0,  # * units.um/units.min,  # when inactive 2-5 um/min, \
        # when active 10-15 um/min (Lammermann, 2008)
        'diffusion': DIFFUSION_RATES,
        'pi': PI,
        'nAvagadro': AVOGADRO / units.mol,  # count / mol, #TODO convert back from ng
        'external_concentration_unit': CONCENTRATION_UNIT,

        # death rates
        'death_apoptosis': 0.5,  # (Naik, 2008)
        'death_time': 4 * 24 * 60 * 60,  # 10 days (*24*60*60 seconds) (Naik, 2008)

        # division
        'divide_prob': 0.6,  # probability of division 24 hr
        'divide_time': 5 * 24 * 60 * 60,  # 5 days (*24*60*60 seconds) (data)

        # transitions
        'internal_tumor_debris_threshold': 415000,  # 415000 This is in counts # (Yang, 2006)

        # membrane equilibrium amounts
        'PDL1p_PDL1_equilibrium': 5e4,
        'PDL1p_MHCI_equilibrium': 5e4,

        # uptake
        'tumor_debris_uptake': 300 / 60,  # number of tumor debris molecules/cell/hr uptaken \
        # conv to seconds (Yang, 2006)
        'tumor_debris_MW': 29000 * units.g / units.mol,  # g/mol (Apetoh, 2007)
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        diffusion_area = self.parameters['diffusion']['tumor_debris'] * self.parameters['timestep'] * units.s
        diffusion_radius = diffusion_area ** 0.5
        sphere_radius = self.parameters['diameter'] / 2 + diffusion_radius
        external_tumor_debris_available_volume = 4 / 3 * self.parameters['pi'] * sphere_radius ** 3
        self.molar_available_tumor_debris = (self.parameters['external_concentration_unit'] * \
                                            external_tumor_debris_available_volume * self.parameters['nAvagadro'] / \
                                            self.parameters['tumor_debris_MW']).to('count').magnitude

    def initial_state(self, config=None):
        return {
            'boundary': {
                'cell_type': 'dendritic'}}

    def ports_schema(self):
        return {
            'globals': {
                'death': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'},  # nice to have. Need mechanism that turns death on. Low probability threshold?
                'divide': {
                    '_default': False,
                    '_updater': 'set'},
                'divide_count': {
                    '_default': 0,
                    '_updater': 'accumulate'}},  # used to count number of divisions over time.
            'internal': {
                'cell_state': {
                    '_default': 'inactive',
                    '_emit': True,
                    '_updater': 'set',
                },  # either 'activate' or 'inactive'
                'tumor_debris': {
                    '_default': 0,
                    '_updater': 'accumulate',
                    '_emit': True},
                'cell_state_count': {
                    '_default': 0,
                    '_updater': 'accumulate'},  # counts how many total cell in a given time. Might not be needed.
                'lymph_node_timer': {
                    # counts how long in lymph node, high value increases chance of migration back to tumor
                    '_default': 0,  # TODO -- does the LN time this, or do the cells?
                    '_emit': True,
                    '_updater': 'accumulate'},
            },
            'boundary': {
                'cell_type': {
                    '_value': 'dendritic',
                    '_emit': True, },
                # Might be needed for neighbors, but really for the experimenters to quantify
                'mass': {
                    '_value': self.parameters['mass']},
                'diameter': {
                    '_default': self.parameters['diameter']},
                'velocity': {
                    '_default': self.parameters['velocity'],
                    '_emit': True},
                'external': {
                    'tumor_debris': {
                        '_default': 0.0,  # TODO: units.ng / units.mL
                        '_emit': True},
                    'lymph_node': {
                        '_default': False}}},  # this is True when in the lymph node, begins counter for how long.
            'neighbors': {  # this is only for presenting in the lymph node, not in the tumor "arena"
                'present': {
                    'PDL1': {
                        '_default': 0,
                        '_updater': 'set'},
                    'MHCI': {
                        '_default': 0,  # high level for activation, should come from environment
                        '_updater': 'set',
                        '_emit': True}},
                'accept': {
                    'PD1': {'_default': 0},
                    'TCR': {
                        '_default': 0,
                        '_emit': True
                    }
                }
            }}

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        external_tumor_debris = states['boundary']['external']['tumor_debris']  # concentration
        internal_tumor_debris = states['internal']['tumor_debris']  # counts

        # determine available tumor debris
        available_tumor_debris_counts = external_tumor_debris * self.molar_available_tumor_debris
        # TODO -- we need to move available tumor debris to internal tumor debris
        # TODO - test out with experiment and also do calculation

        # death by apoptosis
        prob_death = get_probability_timestep(
            self.parameters['death_apoptosis'],
            self.parameters['death_time'],
            timestep)
        if random.uniform(0, 1) < prob_death:
            return {
                'globals': {
                    'death': 'apoptosis'}}

        # division
        if cell_state == 'active':
            prob_divide = get_probability_timestep(
                self.parameters['divide_prob'],
                self.parameters['divide_time'],
                timestep)
            if random.uniform(0, 1) < prob_divide:
                PDL1n_divide_count = 1
                return {
                    'globals': {
                        'divide': True,
                        'PDL1n_divide_count': PDL1n_divide_count
                    }}

        ## Build up an update
        update = {
            'internal': {},
            'boundary': {},
            'neighbors': {
                'present': {}}}

        # state transition
        new_cell_state = cell_state
        if cell_state == 'inactive':
            if internal_tumor_debris >= self.parameters['internal_tumor_debris_threshold']:  # TODO -- the parameter is in ng/mL, while the variable is ints
                new_cell_state = 'active'
                cell_state_count = 1
                update.update({
                    'internal': {
                        'cell_state': new_cell_state,
                        'cell_state_count': cell_state_count}})

        # behavior
        # uptake locally available tumor debris in the environment
        tumor_debris_uptake = min(
            int(self.parameters['tumor_debris_uptake'] * timestep),
            int(available_tumor_debris_counts))  # TODO -- check this
        update['boundary'].update(
            {'exchange': {'tumor_debris': -tumor_debris_uptake}})
        update['internal'].update({
            'tumor_debris': tumor_debris_uptake})

        if new_cell_state == 'active':
            PDL1 = self.parameters['PDL1p_PDL1_equilibrium']
            MHCI = self.parameters['PDL1p_MHCI_equilibrium']
            update['neighbors']['present'].update({
                'PDL1': PDL1,
                'MHCI': MHCI})

        return update


# test functions
def get_timeline(
        total_time=129600,
        number_steps=10):
    """Make a timeline that feeds input to the tumor process"""

    interval = total_time / (number_steps * TIMESTEP)

    timeline = [
        (interval * 0 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 1 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 5e4,
        }),
        (interval * 2 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 3 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 4 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 5 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 6 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 7 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
        (interval * 8 * TIMESTEP, {
            ('boundary', 'external', 'tumor_debris'): 1e7,
            ('neighbors', 'accept', 'TCR'): 0.0,
        }),
    ]
    return timeline


def test_single_d_cell(
        total_time=129600,
        time_step=TIMESTEP,
        out_dir='out'
):
    timeline = get_timeline()

    d_cell_process = DendriticCellProcess({'timestep': time_step})
    timeline_process = TimelineProcess({'timeline': timeline})

    # initial state
    initial_state = d_cell_process.initial_state()
    initial_state['boundary'] = {'external': {'tumor_debris': 0.1}}

    # put in composite
    processes = {
        'd_cell': d_cell_process,
        'timeline': timeline_process
    }
    topology = {
        'd_cell': {port: (port,) for port in d_cell_process.ports_schema().keys()},
        'timeline': {port: (port,) for port in timeline_process.ports()}
    }

    composite = Composite({
        'processes': processes,
        'topology': topology,
        'initial_state': initial_state})

    # run experiment
    sim = Engine(composite=composite)
    sim.update(total_time)

    # get the data
    timeseries = sim.emitter.get_timeseries()

    # plot
    plot_settings = {'remove_zeros': False}
    plot_simulation_output(timeseries, plot_settings, out_dir, 'dendritic_cell_single')


if __name__ == '__main__':
    test_single_d_cell()
