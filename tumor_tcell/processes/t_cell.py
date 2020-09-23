from __future__ import absolute_import, division, print_function

import os
import sys
import copy
import math
import random
import argparse

from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    plot_simulation_output,
    plot_agents_multigen,
)

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'T_cell'
TIMESTEP = 60 # seconds


def get_probability_timestep(probability_parameter, timescale, timestep):
    ''' transition probability as function of time '''
    rate = -math.log(1 - probability_parameter)
    timestep_fraction = timestep / timescale
    return 1 - math.exp(-rate * timestep_fraction)



class TCellProcess(Process):
    """T-cell process with 2 states

    States:
        - PD1p (PD1+)
        - PD1n (PD1-)

    TODOs
        - make this work!
    """

    name = NAME
    defaults = {
        'time_step': TIMESTEP,
        'diameter': 10 * units.um,
        'initial_PD1n': 0.8,
        'transition_PD1n_to_PD1p_10days': 0.95,  # nearly 95% by 10 days (Wherry 2007)
        'cellstate_transition_time': 259200, # extended contact with tumor cells expressing MHCI for more
               # than 3 days (expected to be 100% in tumor/chronic infection after 1 week) (Schietinger 2016)

        # death rates (Petrovas 2007)
        'PDL1_critical_number': 1e4,  # threshold number of PDL1 molecules/cell to trigger death
        'death_PD1p_14hr': 0.7,  # 0.7 / 14 hrs
        'death_PD1n_14hr': 0.2,  # 0.2 / 14 hrs
        'death_PD1p_next_to_PDL1p_14hr': 0.95,  # 0.95 / 14 hrs

        # production rates
        'PD1n_IFNg_production': 1.62e4/3600,  # molecules/cell/second (Bouchnita 2017)
        # TODO @Eran - the other IFNg is in ng/mL how does this production get converted?

        'PD1p_IFNg_production': 0.0,  # molecules/cell/second
        'PD1p_PD1_equilibrium': 5e4,  # equilibrium value of PD1 for PD1p (TODO -- get reference)

        'Ligand_threshold': 1e4,  # molecules/neighbor cell

        # division rate (Petrovas 2007, Vodnala 2019)
        'PD1n_growth_8hr': 0.90,  # 95% division in 28 hours
        'PD1p_growth_8hr': 0.20,  # 20% division in 28 hours

        # migration
        'PD1n_migration': 10.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor': 2.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor_dwell_time': 25.0,  # minutes (Thibaut 2020)
        'PD1p_migration': 5.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor': 1.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor_dwell_time': 10.0,  # minutes (Thibaut 2020)

        # killing  ()
        # TODO - @Eran - This value needs to be multiplied by 100 to deal with timestep usually 0.4 packet/min
        # linear production over 4 hr up to a total of 102+-20 granules # (Betts, 2004), (Zhang, 2006)
        'PD1n_cytotoxic_packets': 40,  # number of packets/min to each contacted tumor cell

        # 1:10 fold reduction of PD1+ T cell cytotoxic production (Zelinskyy, 2005)
        'PD1p_cytotoxic_packets': 4,  # number of packets/min to each contacted tumor cell

        # 4 fold reduction in production in T cells in contact with MHCI- tumor
        # (Bohm, 1998), (Merritt, 2003)
        'MHCIn_reduction_production': 4,

        # settings
        'self_path': tuple(),
    }

    def __init__(self, parameters=None):
        super(TCellProcess, self).__init__(parameters)

        if random.uniform(0, 1) < self.parameters['initial_PD1n']:
            self.initial_state = 'PD1n'
        else:
            self.initial_state = 'PD1p'
        self.self_path = self.parameters['self_path']

    def ports_schema(self):
        return {
            'globals': {
                'death': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'},
                'divide': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'},
                'PD1n_divide_count': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'},
                'PD1p_divide_count': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'}
            },
            'internal': {
                'cell_state': {
                    '_default': self.initial_state,
                    '_emit': True,
                    '_updater': 'set'
                },
                'cell_state_count': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'}
            },
            'boundary': {
                'diameter': {
                    '_default': self.parameters['diameter']
                },
                'IFNg': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },
                'PD1': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'set',
                },  # membrane protein, promotes T-cell death
                'cytotoxic_packets': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # release into the tumor cells
                'MHCI_timer': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # affects transition rate to PD1+
            },
            'neighbors': {
                'PDL1': {
                    '_default': 0,
                    '_emit': True,
                },
                'MHCI': {
                    '_default': 0,
                    '_emit': True,
                },
            }
        }

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        PDL1 = states['neighbors']['PDL1']
        MHCI = states['neighbors']['MHCI']
        MHCI_timer = states['boundary']['MHCI_timer']

        # death
        if cell_state == 'PD1n':
            prob_death = get_probability_timestep(
                self.parameters['death_PD1n_14hr'],
                50400,  # 14 hours (14*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_death:
                # print('DEATH PD1- cell!')
                return {
                    '_delete': {
                        'path': self.self_path
                    },
                    'globals': {
                        'death': 'PD1n_apoptosis',
                    }
                }

        elif cell_state == 'PD1p':
            if PDL1 >= self.parameters['PDL1_critical_number']:
                prob_death = get_probability_timestep(
                    self.parameters['death_PD1p_next_to_PDL1p_14hr'],
                    50400,  # 14 hours (14*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_death:
                    # print('DEATH PD1+ cell with PDL1!')
                    return {
                        '_delete': {
                            'path': self.self_path
                        },
                        'globals': {
                            'death': 'PD1p_PDL1_death',
                        }
                    }

            else:
                prob_death = get_probability_timestep(
                    self.parameters['death_PD1p_14hr'],
                    50400,  # 14 hours (14*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_death:
                    # print('DEATH PD1+ cell without PDL1!')
                    return {
                        '_delete': {
                            'path': self.self_path
                        },
                        'globals': {
                            'death': 'PD1p_apoptosis',
                        }
                    }

        # division
        if cell_state == 'PD1n':
            prob_divide = get_probability_timestep(
                self.parameters['PD1n_growth_8hr'],
                100800,  # 28 hours (28*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_divide:
                # print('DIVIDE PD1- cell!')
                PD1n_divide_count = 1
                return {
                    'globals': {
                        'divide': True,
                        'PD1n_divide_count': PD1n_divide_count
                    }
                }

        elif cell_state == 'PD1p':
            prob_divide = get_probability_timestep(
                self.parameters['PD1p_growth_8hr'],
                100800,  # 28 hours (28*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_divide:
                # print('DIVIDE PD1+ cell!')
                PD1p_divide_count = 1
                return {
                    'globals': {
                        'divide': True,
                        'PD1p_divide_count': PD1p_divide_count
                    }
                }

        update = {}
        # state transition
        new_cell_state = cell_state
        if cell_state == 'PD1n':
            if MHCI_timer > self.parameters['cellstate_transition_time']:
                #print('PD1n become PD1p!')
                new_cell_state = 'PD1p'
                cell_state_count = 1
                if 'internal' not in update:
                    update['internal'] = {}
                update['internal'].update({
                    'cell_state': new_cell_state,
                    'cell_state_count': cell_state_count})

            else:
                cell_state_count = 0
                if MHCI > self.parameters['Ligand_threshold']:
                    if 'internal' not in update:
                        update['internal'] = {}
                    update['internal'].update({
                        'cell_state': new_cell_state,
                        'cell_state_count': cell_state_count})

                    if 'boundary' not in update:
                        update['boundary'] = {}
                    update['boundary'].update({
                        'MHCI_timer': timestep})

        elif cell_state == 'PDL1p':
            cell_state_count = 0

            if 'internal' not in update:
                update['internal'] = {}
            update['internal'].update({
                'cell_state': new_cell_state,
                'cell_state_count': cell_state_count})

        # if cell_state == 'PD1n':
        #     prob_transition = get_probability_timestep(
        #         self.parameters['transition_PD1n_to_PD1p_10days'],
        #         864000,  # 10 days (60*60*24*10 seconds)
        #         timestep)
        #     if random.uniform(0, 1) < prob_transition:
        #         new_cell_state = 'PD1p'
        #         cell_state_count = 1
        #         update.update({
        #             'internal': {
        #                 'cell_state': new_cell_state}})
        #     else:
        #         cell_state_count = 0
        #
        # elif cell_state == 'PD1p':
        #     cell_state_count = 0

        #import ipdb; ipdb.set_trace()

        # behavior
        IFNg = 0
        PD1 = 0
        cytotoxic_packets = 0

        # TODO migration
        #  #also dependent on MHCI+/PDL1+
        #  #50% bound vs. 30% bound in flow cytometry experiment on low vs. high
        #  #Migration rates above in parameters

        # Killing by passing cytotoxic packets to tumor
        if new_cell_state == 'PD1n':

            # TODO - @Eran - how do we make T cells only interact with 1 of neighbors at a time?
            # and how do we reference this cell - see last elif
            # check conditional 4 fold reduction
            if MHCI >= self.parameters['Ligand_threshold']:

                # Max production for both happens for PD1- T cells in contact with MHCI+ tumor
                cytotoxic_packets = self.parameters['PD1n_cytotoxic_packets'] * timestep

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1n_IFNg_production'] * timestep

                if 'boundary' not in update:
                    update['boundary'] = {}
                update['boundary'].update({
                    'IFNg': IFNg,
                    'cytotoxic_packets': cytotoxic_packets
                    })

            elif MHCI < self.parameters['Ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)
                cytotoxic_packets = self.parameters['PD1n_cytotoxic_packets'] / self.parameters['MHCIn_reduction_production'] * timestep

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1n_IFNg_production'] / self.parameters['MHCIn_reduction_production'] * timestep

                if 'boundary' not in update:
                    update['boundary'] = {}
                update['boundary'].update({
                    'IFNg': IFNg,
                    'cytotoxic_packets': cytotoxic_packets
                    })

            # target behavior 3 contacts required for cell death, 1-4 cells killed/day

            # TODO  - elif T cell is not in contact with tumor (no cytotoxic packets)
            #   continue

        elif new_cell_state == 'PD1p':
            PD1 = self.parameters['PD1p_PD1_equilibrium']

            if MHCI >= self.parameters['Ligand_threshold']:

                # Max production for both happens for T cells in contact with MHCI+ tumor
                cytotoxic_packets = self.parameters['PD1p_cytotoxic_packets'] * timestep

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1p_IFNg_production'] * timestep

                if 'boundary' not in update:
                    update['boundary'] = {}
                update['boundary'].update({
                    'IFNg': IFNg,
                    'PD1': PD1,
                    'cytotoxic_packets': cytotoxic_packets
                    })

            elif MHCI < self.parameters['Ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)
                cytotoxic_packets = self.parameters['PD1p_cytotoxic_packets'] / self.parameters['MHCIn_reduction_production'] * timestep

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1p_IFNg_production'] / self.parameters['MHCIn_reduction_production'] * timestep

                if 'boundary' not in update:
                    update['boundary'] = {}
                update['boundary'].update({
                    'IFNg': IFNg,
                    'PD1': PD1,
                    'cytotoxic_packets': cytotoxic_packets
                    })

            # target behavior 3 contacts required for cell death, 1-4 cells killed/day

            # TODO  - elif T cell is not in contact with tumor (no cytotoxic packets)
            #   continue

        return update



def get_timeline(
        total_time=1296000,
        number_steps=10):

    interval = total_time / (number_steps * TIMESTEP)

    timeline = [
        (interval * 0 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (interval * 1 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (interval * 2 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (interval * 3 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (interval * 4 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (interval * 5 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (interval * 6 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (interval * 7 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (interval * 8 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (interval * 9 * TIMESTEP, {}),
    ]
    return timeline


def test_single_t_cell(
    total_time=43200,
    timeline=None,
    out_dir='out'):

    t_cell_process = TCellProcess({})
    if timeline is not None:
        settings = {
            'timeline': {
                'timeline': timeline}}
    else:
        settings = {'total_time': total_time}

    # run experiment
    timeseries = simulate_process_in_experiment(t_cell_process, settings)

    # plot
    plot_settings = {
        'remove_zeros': False
    }
    plot_simulation_output(timeseries, plot_settings, out_dir, NAME + '_single')


def test_batch_t_cell(
    total_time=43200,
    batch_size=2,
    timeline=None,
    out_dir='out'):

    combined_raw_data = {}
    for single_idx in range(batch_size):
        t_cell_process = TCellProcess({})
        if timeline is not None:
            sim_settings = {
                'timeline': {
                    'timeline': timeline},
                'return_raw_data': True}
        else:
            sim_settings = {
                'total_time': total_time,
                'return_raw_data': True}
        # run experiment
        raw_data = simulate_process_in_experiment(t_cell_process, sim_settings)
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

    parser = argparse.ArgumentParser(description='ODE expression')
    parser.add_argument('--single', '-s', action='store_true', default=False)
    parser.add_argument('--batch', '-b', action='store_true', default=False)
    parser.add_argument('--timeline', '-t', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    total_time = 43200
    if args.timeline or no_args:
        timeline = get_timeline()
        test_single_t_cell(
            timeline=timeline,
            out_dir=out_dir)

    if args.single:
        test_single_t_cell(
            total_time=total_time,
            out_dir=out_dir)

    if args.batch:
        timeline = get_timeline()
        test_batch_t_cell(
            batch_size=10,
            # total_time=300,
            timeline=timeline,
            out_dir=out_dir)