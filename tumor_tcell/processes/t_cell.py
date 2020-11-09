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
from vivarium.core.composition import simulate_process_in_experiment
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium.plots.simulation_output import plot_simulation_output

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'T_cell'
TIMESTEP = 60  # seconds
CONCENTRATION_UNIT = 1  # TODO (ERAN) set value -- units.ng / units.mL


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

    Expected Behavior
        - approx 3 contacts required for cell death, 1-4 cells killed/day
    """

    name = NAME
    defaults = {
        'time_step': TIMESTEP,
        'diameter': 10 * units.um,  # 0.01 * units.mm,
        'initial_PD1n': 0.8,
        'transition_PD1n_to_PD1p_10days': 0.95,  # nearly 95% by 10 days (Wherry 2007)
        'cellstate_transition_time': 259200,  # extended contact with tumor cells expressing MHCI for more
               # than 3 days (expected to be 100% in tumor/chronic infection after 1 week) (Schietinger 2016)

        #Time before TCR downregulation
        'activation_time': 21600,  # activation enables 6 hours of activation and production of cytokines
        # before enters refractory state due to downregulation of TCR (Salerno, 2017), (Gallegos, 2016)
        'activation_refractory_time': 86400,  # refractory period of 18 hours (plus original 6 h activation time)
        # after activation period with limited ability to produce cytokines and cytotoxic packets and to
        # interact with MHCI (Salerno, 2017), (Gallegos, 2016)
        'TCR_downregulated': 0, #reduce TCR to 0 if activated more than 6 hours for another 18 hours
        'TCR_upregulated': 50000,  # reduce TCR to 0 if activated more than 6 hours for another 18 hours

        # death rates (Petrovas 2007)
        'PDL1_critical_number': 1e4,  # threshold number of PDL1 molecules/cell to trigger death
        'death_PD1p_14hr': 0.7,  # 0.7 / 14 hrs
        'death_PD1n_14hr': 0.2,  # 0.2 / 14 hrs
        'death_PD1p_next_to_PDL1p_14hr': 0.95,  # 0.95 / 14 hrs

        # production rates
        'PD1n_IFNg_production': 1.62e4/3600 * CONCENTRATION_UNIT,  # molecules/cell/second (Bouchnita 2017)
        # TODO @Eran - the other IFNg is in ng/mL how does this production get converted?

        'PD1p_IFNg_production': 0.0 * CONCENTRATION_UNIT,  # molecules/cell/second
        'PD1p_PD1_equilibrium': 5e4,  # equilibrium value of PD1 for PD1p (TODO -- get reference)

        'ligand_threshold': 1e4,  # molecules/neighbor cell

        # division rate (Petrovas 2007, Vodnala 2019)
        'PD1n_growth_8hr': 0.80,  # 90% division in 28 hours
        'PD1p_growth_8hr': 0.20,  # 20% division in 28 hours

        # migration
        'PD1n_migration': 10.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor': 2.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor_dwell_time': 25.0,  # minutes (Thibaut 2020)
        'PD1p_migration': 5.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor': 1.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor_dwell_time': 10.0,  # minutes (Thibaut 2020)

        # killing
            # These values need to be multiplied by 100 to deal with timestep usually 0.4 packet/min
        # linear production over 4-6 hr up to a total of 102+-20 granules # (Betts, 2004), (Zhang, 2006)
        'cytotoxic_packet_production': 40/60,  # number of packets/min produced in T cells ##converted to packets/seconds
        'PD1n_cytotoxic_packets_max':10000, #max number able to produce

        # 1:10 fold reduction of PD1+ T cell cytotoxic production (Zelinskyy, 2005)
        'PD1p_cytotoxic_packets_max': 1000,  # max number able to produce

        # 4 fold reduction in production in T cells in contact with MHCI- tumor
        # (Bohm, 1998), (Merritt, 2003)
        'MHCIn_reduction_production': 4,

        #Cytotoxic packet transfer rate for a minute timeperiod
        'cytotoxic_transfer_rate':400, #number of packets/min that can be transferred to tumor cells

        # settings
        'self_path': tuple(),
    }

    def __init__(self, parameters=None):
        super(TCellProcess, self).__init__(parameters)
        self.self_path = self.parameters['self_path']

    def initial_state(self, config=None):
        if random.uniform(0, 1) < self.parameters['initial_PD1n']:
            initial_state = 'PD1n'
        else:
            initial_state = 'PD1p'

        return {
            'internal': {
                'cell_state': initial_state,
            },
            'boundary': {
                'diameter': self.parameters['diameter']
            },
            'neighbors': {
                'present': {
                    'TCR': 50000,
                }
            }
        }

    def ports_schema(self):
        initial_cell_state = 'PD1n' if random.uniform(0, 1) < self.parameters['initial_PD1n'] else 'PD1p'
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
                    '_divider': 'zero',
                    '_updater': 'accumulate'},
                'PD1p_divide_count': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'}
            },
            'internal': {
                'cell_state': {
                    '_default': initial_cell_state,
                    '_emit': True,
                    '_updater': 'set'
                },
                'cell_state_count': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'
                },
                'total_cytotoxic_packets': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate'
                },
                'TCR_timer': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # affects TCR expression
            },
            'boundary': {
                'cell_type': {
                    '_value': 't-cell'
                },
                'diameter': {
                    '_default': self.parameters['diameter'],
                    '_divider': 'set',
                },
                'external': {
                    'IFNg': {
                        '_default': 0.0 * CONCENTRATION_UNIT,
                        '_emit': True,
                        '_updater': 'accumulate',
                    }},
                'MHCI_timer': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # affects transition rate to PD1+
            },
            'neighbors': {
                'present': {
                    'PD1': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'set',
                    }, # membrane protein, promotes T-cell death
                    'TCR': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'set',
                    } # level of TCR that interacts with MHCI on tumor
                },
                'accept': {
                    'PDL1': {
                        '_default': 0,
                        '_updater': 'set',
                        '_emit': True,
                    },
                    'MHCI': {
                        '_default': 0,
                        '_updater': 'set',
                        '_emit': True,
                    }
                },
                'transfer': {
                    'cytotoxic_packets': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'accumulate',
                        '_divider': 'split',
                    }  # release into the tumor cells
                }
            }
        }

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        PDL1 = states['neighbors']['accept']['PDL1']
        MHCI = states['neighbors']['accept']['MHCI']
        MHCI_timer = states['boundary']['MHCI_timer']
        TCR_timer = states['internal']['TCR_timer']
        TCR = states['neighbors']['present']['TCR']

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
                        'path': self.self_path},
                    'globals': {
                        'death': 'PD1n_apoptosis'}}

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
                            'path': self.self_path},
                        'globals': {
                            'death': 'PD1p_PDL1_death'}}

            else:
                prob_death = get_probability_timestep(
                    self.parameters['death_PD1p_14hr'],
                    50400,  # 14 hours (14*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_death:
                    # print('DEATH PD1+ cell without PDL1!')
                    return {
                        '_delete': {
                            'path': self.self_path},
                        'globals': {
                            'death': 'PD1p_apoptosis'}}

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
                        'PD1n_divide_count': PD1n_divide_count}}

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
                        'PD1p_divide_count': PD1p_divide_count}}


        ## Build up an update
        update = {
            'internal': {},
            'boundary': {},
            'neighbors': {'present': {}, 'accept': {}, 'transfer': {}}}

        #TCR downregulation after 6 hours of activation
        if TCR_timer > self.parameters['activation_time']:
            TCR = self.parameters['TCR_downregulated']
            update['neighbors']['present'].update({
                'TCR': TCR})

        #update the TCR_timer for each timestep in contact with MHCI or during refractory period regardless
        if MHCI > self.parameters['ligand_threshold']:
            update['internal'].update({
                'TCR_timer': timestep})
        elif TCR_timer > self.parameters['activation_time']:
            update['internal'].update({
                'TCR_timer': timestep})

        #After the refractory period is over then
        if TCR_timer > self.parameters['activation_refractory_time']:
            TCR = self.parameters['TCR_upregulated']
            update['internal'].update({
                'TCR_timer': -self.parameters['activation_refractory_time']})
            update['neighbors']['present'].update({
                'TCR': TCR})

        # state transition
        new_cell_state = cell_state
        if cell_state == 'PD1n':
            if MHCI_timer > self.parameters['cellstate_transition_time']:
                #print('PD1n become PD1p!')
                new_cell_state = 'PD1p'
                cell_state_count = 1
                update['internal'].update({
                    'cell_state': new_cell_state,
                    'cell_state_count': cell_state_count})

            else:
                cell_state_count = 0
                if MHCI > self.parameters['ligand_threshold']:
                    update['internal'].update({
                        'cell_state': new_cell_state,
                        'cell_state_count': cell_state_count})
                    update['boundary'].update({
                        'MHCI_timer': timestep})

        elif cell_state == 'PD1p':
            cell_state_count = 0
            update['internal'].update({
                'cell_state': new_cell_state,
                'cell_state_count': cell_state_count})


        # behavior
        # TODO migration
        #  #also dependent on MHCI+/PDL1+
        #  #50% bound vs. 30% bound in flow cytometry experiment on low vs. high
        #  #Migration rates above in parameters

        # Killing by passing cytotoxic packets to tumor
        if new_cell_state == 'PD1n':

            # and how do we reference this cell - see last elif
            # check conditional 4 fold reduction
            if MHCI >= self.parameters['ligand_threshold'] and TCR >= self.parameters['ligand_threshold']:

                # Max production for both happens for PD1- T cells in contact with MHCI+ tumor
                #Can max the amount of cytotoxic packets released
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1n_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1n_IFNg_production'] * timestep
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

            elif MHCI > 0 and TCR >= self.parameters['ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)

                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1n_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] / self.parameters[
                        'MHCIn_reduction_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1n_IFNg_production'] / self.parameters['MHCIn_reduction_production'] * timestep
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

            elif MHCI == 0:

                IFNg = 0
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

        elif new_cell_state == 'PD1p':
            PD1 = self.parameters['PD1p_PD1_equilibrium']

            if MHCI >= self.parameters['ligand_threshold'] and TCR >= self.parameters['ligand_threshold']:

                # Max production for both happens for T cells in contact with MHCI+ tumor
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1p_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                update['neighbors']['present'].update({
                    'PD1': PD1})

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1p_IFNg_production'] * timestep
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

            elif MHCI > 0 and TCR >= self.parameters['ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1p_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] / \
                                        self.parameters['MHCIn_reduction_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                update['neighbors']['present'].update({
                    'PD1': PD1})

                # produce IFNg  # rates are determined above
                IFNg = self.parameters['PD1p_IFNg_production'] / \
                       self.parameters['MHCIn_reduction_production'] * timestep
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

            elif MHCI == 0:
                IFNg = 0
                update['boundary'].update({
                    'external': {'IFNg': IFNg}})

        #keep cytotoxic transfer to max limit between tumor and T cells and remove from T cell total when transfer
        #max rate is 400 packets/minute

        max_cytotoxic_transfer = self.parameters['cytotoxic_transfer_rate'] - states['neighbors']['transfer']['cytotoxic_packets']
        cytotoxic_transfer = min(max_cytotoxic_transfer,states['internal']['total_cytotoxic_packets'])
        update['internal'].update({
            'cytotoxic_packets': -cytotoxic_transfer})
        update['neighbors']['transfer'].update({
            'cytotoxic_packets': cytotoxic_transfer})

        return update



def get_timeline(
        total_time=200000,
        number_steps=10):

    interval = total_time / (number_steps * TIMESTEP)

    timeline = [
        (interval * 0 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 0.0,
            ('neighbors', 'accept', 'MHCI'): 0.0,
        }),
        (interval * 1 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 5e5,
            ('neighbors', 'accept', 'MHCI'): 5e5,
        }),
        (interval * 2 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 0.0,
            ('neighbors', 'accept', 'MHCI'): 0.0,
        }),
        (interval * 3 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 0.0,
            ('neighbors', 'accept', 'MHCI'): 0.0,
        }),
        (interval * 4 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 0.0,
            ('neighbors', 'accept', 'MHCI'): 0.0,
        }),
        (interval * 5 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 5e5,
            ('neighbors', 'accept', 'MHCI'): 5e5,
        }),
        (interval * 6 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 5e5,
            ('neighbors', 'accept', 'MHCI'): 5e5,
        }),
        (interval * 7 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 0.0,
            ('neighbors', 'accept', 'MHCI'): 0.0,
        }),
        (interval * 8 * TIMESTEP, {
            ('neighbors', 'accept', 'PDL1'): 5e5,
            ('neighbors', 'accept', 'MHCI'): 5e5,
        }),
        (interval * 9 * TIMESTEP, {}),
    ]
    return timeline


def test_single_t_cell(
    total_time=43200,
    time_step=TIMESTEP,
    timeline=None,
    out_dir='out'):

    t_cell_process = TCellProcess({})
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
    settings['initial_state'] = t_cell_process.initial_state()

    # run experiment
    timeseries = simulate_process_in_experiment(t_cell_process, settings)

    # plot
    plot_settings = {
        'remove_zeros': False
    }
    plot_simulation_output(timeseries, plot_settings, out_dir, NAME + '_single')


def test_batch_t_cell(
    total_time=43200,
    time_step=TIMESTEP,
    batch_size=2,
    timeline=None,
    out_dir='out'):

    combined_raw_data = {}
    for single_idx in range(batch_size):
        t_cell_process = TCellProcess({})
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
        sim_settings['initial_state'] = t_cell_process.initial_state()

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
