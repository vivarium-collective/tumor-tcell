from __future__ import absolute_import, division, print_function

import os
import sys
import math
import random
import argparse

from vivarium.library.units import units
from vivarium.core.process import Process
from vivarium.core.composition import simulate_process
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium.plots.simulation_output import plot_simulation_output

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'Tumor'
TIMESTEP = 60
CONCENTRATION_UNIT = 1  # TODO (ERAN) set value -- units.ng / units.mL


def get_probability_timestep(probability_parameter, timescale, timestep):
    ''' transition probability as function of time '''
    rate = -math.log(1 - probability_parameter)
    timestep_fraction = timestep / timescale
    return 1 - math.exp(-rate * timestep_fraction)


class TumorProcess(Process):
    """Tumor process with 2 states

    States:
        - PDL1p (PDL1+, MHCI+)
        - PDL1n (PDL1-, MHCI-)

    """

    name = NAME
    defaults = {
        'time_step': TIMESTEP,
        'diameter': 15 * units.um, #  0.01 * units.mm,
        'mass': 8 * units.ng,
        'initial_PDL1n': 0.9, #most start out this way based on data

        # death rates
        'death_apoptosis': 0.5,  # negligible compared to growth/killing 0.95 by 5 day (Gong, 2017)
        'cytotoxic_packet_threshold': 128 * 100,  # need at least 128 packets for death (multiply by 10 for T cells)

        # division rate
        'PDL1n_growth': 0.6,  # probability of division 24 hr (Eden, 2011)
        #'PDL1p_growth': 0,  # Cells arrested - do not divide (data, Thibaut 2020, Hoekstra 2020)

        # migration
        #'tumor_migration': 0.25,  # um/minute (Weigelin 2012) #set to 0 for now because smaller than T cells

        #IFNg Internalization
        'Max_IFNg_internalization': 21/60, #number of IFNg 1250 molecules/cell/hr degraded conv to seconds
        #volume to convert counts to available IFNg molecules able to be internalized based on the diffusion
        # coefficient and timestep of 60s
        'external_IFNg_available_volume': 8.24*10 ** -8, #* units.mL, # in mL 12 um +diameter of 15 um = 4/3*pi*(27 um)^3
        #TODO - make this more general from timestep/diameter
        #TODO - Use a global IFNg MW (also use in local fields)
        #TODO - synchronize expected concentration units with local fields
        'Avagadro_num': 6.022*10 ** 14, #* units.count / units.nmol, #convert back from ng
        'IFNg_MW': 17000, #* units.ng / units.nmol,
        'IFNg_threshold': 15000, #calculated from home data of incubating 1 ng/mL for 20 mL and 20x10^6 cells and half-life
        'reduction_IFNg_internalization': 2, #based on data from (Ersvaer, 2007) & (Darzi, 2017)

        # membrane equilibrium amounts
        'PDL1p_PDL1_equilibrium': 5e4, #TODO ref
        'PDL1p_MHCI_equilibrium': 5e4, #TODO ref

        # settings
        'self_path': tuple(),
    }

    def __init__(self, parameters=None):
        super(TumorProcess, self).__init__(parameters)
        self.self_path = self.parameters['self_path']

    def initial_state(self, config=None):
        if random.uniform(0, 1) < self.parameters['initial_PDL1n']:
            initial_state = 'PDL1n'
        else:
            initial_state = 'PDL1p'

        return {
            'internal': {
                'cell_state': initial_state,
            },
            'boundary': {
                'diameter': self.parameters['diameter']
            },
        }

    def ports_schema(self):
        initial_cell_state = 'PDL1n' if random.uniform(0, 1) < self.parameters['initial_PDL1n'] else 'PDL1p'
        return {
            'globals': {
                'death': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'},
                'divide': {
                    '_default': False,
                    #'_emit': False, #true for monitoring behavior in process
                    '_updater': 'set'},
                'PDL1n_divide_count': {
                    '_default': 0,
                    #'_emit': False, #true for monitoring behavior in process
                    '_divider': 'zero',
                    '_updater': 'accumulate'}
            },
            'internal': {
                'cell_state': {
                    '_default': initial_cell_state,
                    '_emit': True,
                    '_updater': 'set'
                },
                'IFNg': {
                    '_default': 0,
                    '_emit': True,
                    '_divider': 'split',
                    '_updater': 'accumulate'
                },
                'cell_state_count': {
                    '_default': 0,
                    #'_emit': False, #true for monitoring behavior in process
                    '_updater': 'accumulate'}
            },
            'boundary': {
                'cell_type': {
                    '_value': 'tumor'
                },
                'mass': {
                    '_value': self.parameters['mass']
                },
                'diameter': {
                    '_default': self.parameters['diameter']
                },
                'velocity': {
                    '_default': 0.0 * units.um / units.s,
                    '_updater': 'set',
                },
                'external': {
                    'IFNg': {
                        '_default': 0.0 * CONCENTRATION_UNIT,
                        '_emit': True,
                    }},  # cytokine changes tumor phenotype to MHCI+ and PDL1+
                'exchange': {
                    'IFNg': {
                        '_default': 0,  # counts
                        # '_emit': False, #true for monitoring behavior in process
                        '_updater': 'accumulate',
                    }},
            },
            'neighbors': {
                'present': {
                    'PDL1': {
                        '_default': 0,
                        #'_emit': False, #true for monitoring behavior in process
                        '_updater': 'set',
                    },  # membrane protein, promotes T cell exhuastion and deactivation with PD1
                    'MHCI': {
                        '_default': 1000,
                        #'_emit': False, #true for monitoring behavior in process
                        '_updater': 'set',
                    }  # membrane protein, promotes Tumor death and T cell activation with TCR
                },
                'accept': {
                    'PD1': {
                        '_default': 0,
                        #'_emit': False, #true for monitoring behavior in process
                        '_updater': 'set',
                    },
                    'TCR': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'set',
                    }
                },
                'receive': {
                    'cytotoxic_packets': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'accumulate',
                        '_divider': 'split',
                    } #from T cells
                }
            }
        }

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        cytotoxic_packets = states['neighbors']['receive']['cytotoxic_packets']
        external_IFNg = states['boundary']['external']['IFNg'] #concentration
        internal_IFNg = states['internal']['IFNg'] #counts
        #IFNg_timer = states['boundary']['IFNg_timer']


        #determine available IFNg
        available_IFNg = external_IFNg * self.parameters['external_IFNg_available_volume'] \
                         * self.parameters['Avagadro_num'] / self.parameters['IFNg_MW']

        # death by apoptosis
        prob_death = get_probability_timestep(
            self.parameters['death_apoptosis'],
            432000,  # 5 days (5*24*60*60 seconds)
            timestep)
        if random.uniform(0, 1) < prob_death:
            #print('Apoptosis!')
            return {
                'globals': {
                    'death': 'apoptosis'}}


        if cytotoxic_packets >= self.parameters['cytotoxic_packet_threshold']:
            #print('Tcell_death!')
            return {
                'globals': {
                    'death': 'Tcell_death'}}

        # division
        if cell_state == 'PDL1n':
            prob_divide = get_probability_timestep(
                self.parameters['PDL1n_growth'],
                86400,  # 24 hours (24*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_divide:
                #print('PDL1n DIVIDE!')
                PDL1n_divide_count = 1
                return {
                    'globals': {
                        'divide': True,
                        'PDL1n_divide_count': PDL1n_divide_count
                    }
                }
        elif cell_state == 'PDL1p':
            pass


        ## Build up an update
        update = {'internal': {},
                  'boundary': {},
                  'neighbors': {'present': {}, 'accept': {}, 'receive': {}}}

        # state transition
        new_cell_state = cell_state
        if cell_state == 'PDL1n':
            if internal_IFNg >= self.parameters['IFNg_threshold']:
                new_cell_state = 'PDL1p'
                cell_state_count = 1
                update.update({
                    'internal': {
                        'cell_state': new_cell_state,
                        'cell_state_count': cell_state_count}})

        elif cell_state == 'PDL1p':
            cell_state_count = 0

        # behavior
        MHCI = 1000
        PDL1 = 0

        if new_cell_state == 'PDL1p':
            PDL1 = self.parameters['PDL1p_PDL1_equilibrium']
            MHCI = self.parameters['PDL1p_MHCI_equilibrium']

            update['neighbors']['present'].update({
                'PDL1': PDL1,
                'MHCI': MHCI})

            # degrade locally available IFNg in the environment
            IFNg_degrade = min(int(self.parameters['Max_IFNg_internalization'] \
                                   / self.parameters['reduction_IFNg_internalization'] * timestep), int(available_IFNg))

            update['boundary'].update({
                'exchange': {'IFNg': -IFNg_degrade}})
            update['internal'].update({
                'IFNg': IFNg_degrade})

        elif new_cell_state == 'PDL1n':
            #degrade locally available IFNg in the environment
            IFNg_degrade = min(int(self.parameters['Max_IFNg_internalization'] * timestep), int(available_IFNg))

            update['boundary'].update({
                'exchange': {'IFNg': -IFNg_degrade}})
            update['internal'].update({
                'IFNg': IFNg_degrade})


        return update



def get_timeline(
        total_time=129600,
        number_steps=10):

    interval = total_time/(number_steps*TIMESTEP)

    timeline = [
        (interval * 0 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 0.0,
            ('boundary', 'external', 'IFNg'): 0.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 0.0,
        }),
        (interval * 1 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 1000.0,
            ('boundary', 'external', 'IFNg'): 1.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 2 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 2000.0,
            ('boundary', 'external', 'IFNg'): 2.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 0.0,
        }),
        (interval * 3 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 3000.0,
            ('boundary', 'external', 'IFNg'): 3.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 4 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 4000.0,
            ('boundary', 'external', 'IFNg'): 4.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 5 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 7000.0,
            ('boundary', 'external', 'IFNg'): 3.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 6 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 10000.0,
            ('boundary', 'external', 'IFNg'): 2.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 7 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 15000.0,
            ('boundary', 'external', 'IFNg'): 2.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 8 * TIMESTEP, {
            ('neighbors', 'receive', 'cytotoxic_packets'): 16000.0,
            ('boundary', 'external', 'IFNg'): 2.0 * CONCENTRATION_UNIT,
            ('neighbors', 'accept', 'PD1'): 5e4,
        }),
        (interval * 9 * TIMESTEP, {}),
    ]
    return timeline


def test_single_Tumor(
        total_time=43200,
        time_step=TIMESTEP,
        timeline=None,
        out_dir='out'):

    Tumor_process = TumorProcess({})

    if timeline is not None:
        settings = {
            'timeline': {
                'timeline': timeline,
                'time_step': time_step}}
    else:
        settings = {
            'total_time': total_time,
            'time_step': time_step}

    # get initial state
    settings['initial_state'] = Tumor_process.initial_state()

    # run experiment
    timeseries = simulate_process(Tumor_process, settings)

    # plot
    plot_settings = {'remove_zeros': False}
    plot_simulation_output(timeseries, plot_settings, out_dir, NAME + '_single')

def test_batch_tumor(
    total_time=43200,
    time_step=TIMESTEP,
    batch_size=2,
    timeline=None,
    out_dir='out'):

    combined_raw_data = {}
    for single_idx in range(batch_size):
        Tumor_process = TumorProcess({})
        if timeline is not None:
            sim_settings = {
                'timeline': {
                    'timeline': timeline,
                    'time_step': time_step},
                'return_raw_data': True}
        else:
            sim_settings = {
                'total_time': total_time,
                'time_step': time_step,
                'return_raw_data': True}

        # get initial state
        sim_settings['initial_state'] = Tumor_process.initial_state()

        # run experiment
        raw_data = simulate_process(Tumor_process, sim_settings)
        for time, time_data in raw_data.items():
            if time not in combined_raw_data:
                combined_raw_data[time] = {'agents': {}}
            combined_raw_data[time]['agents'][single_idx] = time_data

    plot_settings = {
        'agents_key': 'agents'}
    plot_agents_multigen(combined_raw_data, plot_settings, out_dir, NAME + '_batch')


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='tumor cells')
    parser.add_argument('--single', '-s', action='store_true', default=False)
    parser.add_argument('--batch', '-b', action='store_true', default=False)
    parser.add_argument('--timeline', '-t', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    total_time = 43200
    if args.single:
        test_single_Tumor(
            total_time=total_time,
            out_dir=out_dir)

    if args.timeline or no_args:
        timeline = get_timeline()
        test_single_Tumor(
            timeline=timeline,
            out_dir=out_dir)

    if args.batch:
        timeline = get_timeline()
        test_batch_tumor(
            batch_size=10,
            # total_time=300,
            timeline=timeline,
            out_dir=out_dir)