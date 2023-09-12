"""
==============
T Cell Process
==============

The T cell process is focused on two states of a T cell: PD1- with increased
secretion of immune molecules (IFNg and cytotoxic packets) and PD1+ with
decreased secretion of immune molecules (IFNg and cytotoxic packets). These
immune molecules have impact of the state and death of tumor cells. Its
transition from the PD1- state is dependent on the length of time it is engaged
with tumor cells.

This process can be run on its own from the command line. There are three simulation
options: "single", "batch" and "timeline". Single simulates the process on its own
one time, batch simulates the process multiple times to demonstrate the stochasticity,
and timeline simulate the process one time with pre-set perturbations that mimic
interactions with other cells.

    $ python tumor_tcell/processes/t_cell.py [--single, -s] [--batch, -b] [--timeline, -t]


"""

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


NAME = 'T_cell'
TIMESTEP = 60  # seconds
CONCENTRATION_UNIT = 1  # TODO: units.ng / units.mL

def lymph_node_division(mother_value, **args):
    """Sets mother cell to true"""
    if mother_value:
        return [True, False]
    else:
        return [False, False]

def assymetric_division(mother_value, **args):
    """One daughter value's count increases"""
    return [mother_value+1, mother_value]

def set_velocity_default(mother_value, **args):
    """Resets daughter velocities"""
    return [10.0 * units.um/units.min, 10.0 * units.um/units.min]

def get_probability_timestep(probability_parameter, timescale, timestep):
    """ transition probability as function of time """
    rate = -math.log(1 - probability_parameter)
    timestep_fraction = timestep / timescale
    return 1 - math.exp(-rate * timestep_fraction)

def probability_of_occurrence_within_interval(interval_duration, expected_time):
    """
    Compute the probability that an event will occur at least once
    within a given time interval.

    This assumes the event follows a Poisson process, where the event
    is expected to occur once every `expected_time` hours.

    Args:
        interval_duration (float): The duration of the time interval.
        expected_time (float): The expected time between occurrences.

    Returns:
        float: The probability of the event occurring at least once
               within the time interval.
    """
    lambda_ = interval_duration / expected_time
    P_0 = math.exp(-lambda_)
    P_at_least_one = 1 - P_0
    return P_at_least_one


# def get_probability_timestep(probability_parameter, interval_duration, expected_time):
#     """
#     Calculate the probability of an event occurring within a given
#     time interval, given the average time between occurrences and
#     the probability of an event occurring within that average time.
#
#     Args:
#         interval_duration (float): The duration of the time interval
#                                    in hours.
#         expected_time (float): The average time between occurrences
#                                of the event in hours.
#         probability_parameter (float): The probability of the event occurring
#                              within the expected time.
#
#     Returns:
#         float: The probability of the event occurring within the
#                interval_duration.
#     """
#     # Compute the average rate of occurrence
#     lambda_ = -math.log(1 - probability_parameter) / expected_time
#
#     # Calculate the average rate of occurrence within the desired interval
#     lambda_interval = lambda_ * interval_duration
#
#     # Compute the probability of no occurrence
#     P_0 = math.exp(-lambda_interval)
#
#     # The probability of at least one occurrence is 1 - P(0)
#     P_at_least_one = 1 - P_0
#
#     return P_at_least_one


class TCellProcess(Process):
    """TCellProcess

    References:
        * Salerno, F., Paolini, N. A., Stark, R., von Lindern, M., & Wolkers, M. C. (2017). Distinct PKC-mediated posttranscriptional events set cytokine production kinetics in CD8+ T cells. Proceedings of the National Academy of Sciences, 114(36), 9677-9682.
    """
    name = NAME
    defaults = {
        'time_step': TIMESTEP,
        'diameter': 7.5 * units.um,
        'mass': 2 * units.ng,
        'initial_PD1n': 0.8,
        'refractory_count_threshold': 3,  # assuming that cells have been stimulated twice already in culture
        'external_molecules': ['IFNg', 'tumor_debris'],
        # and need 5 stimulations to become exhausted (Zhao, 2020)

        #Time before TCR downregulation
        'activation_time': 21600,  # activation enables 6 hours of activation and production of cytokines
        # before enters refractory state due to downregulation of TCR (Salerno, 2017), (Gallegos, 2016)
        'activation_refractory_time': 43200,  # refractory period of 18 hours (plus original 6 h activation time)
        # after activation period with limited ability to produce cytokines and cytotoxic packets and to
        # interact with MHCI (Salerno, 2017), (Gallegos, 2016)
        'TCR_downregulated': 0, # reduce TCR to 0 if activated more than 6 hours for another 18 hours
        'TCR_upregulated': 50000,  # restore TCR after refractory period

        # death rates
        'PDL1_critical_number': 1e4,  # threshold number of PDL1 molecules/cell to induce behavior
        'death_PD1p_14hr': 0.35,  # 0.7 / 14 hrs (Petrovas 2007)
        'death_PD1n_14hr': 0.1,  # 0.2 / 14 hrs (Petrovas 2007)
        'death_PD1p_next_to_PDL1p_14hr': 0.475,  # 0.95 / 14 hrs (Dong 2002, Tang 2015)

        # production rates
        'PD1n_IFNg_production': 1.62e4/3600,  # molecule counts/cell/second (Bouchnita 2017)
        'PD1p_IFNg_production': 1.62e3/3600,  # molecule counts/cell/second (Zelinskyy, 2005)
        'PD1p_PD1_equilibrium': 5e4,  # equilibrium value of PD1 for PD1p

        'ligand_threshold': 1e4,  # molecules/neighbor cell needed to recognize ligands on positive cells

        # division rate
        'PD1n_growth_28hr': 0.90,  # 90% division in 28 hours (Petrovas 2007, Vodnala 2019)
        'PD1p_growth_28hr': 0.20,  # 20% division in 28 hours (Petrovas 2007, Vodnala 2019)
        'PD1n_divide_threshold': 5,  # counts for triggering division (Zhao, 2020)
        'LymphNode_delay_growth': 32400, # 14400 divide approximately every 4 hours, or 5-6 times in 24 hours. \
            # 4*60*60=14400, (Mempel, 2004);(Bousso, 2008)


        # migration
        'PD1n_migration': 10.0 * units.um/units.min,  # um/minute (Boissonnas 2007)
        'migration_MHCIp_tumor_dwell_velocity': 0.0 * units.um/units.min, #(Thibaut 2020)
        'PD1n_migration_MHCIp_tumor_dwell_time': 25.0*60,  # minutes converted to seconds (Thibaut 2020)
        'PD1p_migration': 5.0 * units.um/units.min,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor_dwell_time': 10.0*60,  # minutes converted to seconds (Thibaut 2020)
        'PD1n_migration_refractory_time': 35.0 * 60,  # 25 minutes of refractory where not interact with tumor (Thibaut 2020)
        'PD1p_migration_refractory_time': 20.0 * 60,  # 10 minutes of refractory where not interact with tumor (Thibaut 2020)

        # killing
        # These values need to be multiplied by 100 to deal with timestep usually 0.4 packet/min
        # linear production over 4-6 hr up to a total of 102+-20 granules (Betts, 2004), (Zhang, 2006)
        'cytotoxic_packet_production': 40/60,  # number of packets/min produced in T cells ##converted to packets/seconds
        'PD1n_cytotoxic_packets_max': 10000,  # max number able to produce
        # 1:10 fold reduction of PD1+ T cell cytotoxic production (Zelinskyy, 2005)
        'PD1p_cytotoxic_packets_max': 1000,  # max number able to produce
        # 4 fold reduction in production in T cells in contact with MHCI- tumor (Bohm, 1998), (Merritt, 2003)
        'MHCIn_reduction_production': 400,
        # Cytotoxic packet transfer rate for a minute time period
        'cytotoxic_transfer_rate': 400,  #number of packets/min that can be transferred to tumor cells (Betts, 2004), (Zhang, 2006)
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

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
                    '_updater': 'set'},
                'PD1n_divide_count': {
                    '_default': 0,
                    '_divider': {
                        'divider': assymetric_division,
                        'config': {}
                    },
                    '_updater': 'accumulate'},
                'PD1p_divide_count': {
                    '_default': 0,
                    '_updater': 'accumulate'},
                'LN_no_migration': {   # TODO (ERAN) -- can this be removed??
                    '_default': False,
                    '_divider': {
                        'divider': lymph_node_division,
                        'config': {}}},
            },
            'internal': {
                'cell_state': {
                    '_default': initial_cell_state,
                    '_emit': True,
                    '_updater': 'set'
                },
                'cell_state_count': {
                    '_default': 0,
                    '_updater': 'accumulate'
                },
                'refractory_count': {
                    '_default': 0,
                    '_updater': 'accumulate'
                },
                'total_cytotoxic_packets': {
                    '_default': 0,
                    '_updater': 'accumulate',
                    '_divider': 'split',
                },
                'TCR_timer': {
                    '_default': 0,
                    '_updater': 'accumulate',
                },  # affects TCR expression
                'velocity_timer': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # affects dwell time at tumor
            },
            'boundary': {
                'cell_type': {
                    '_value': 't-cell'
                },
                'mass': {
                    '_value': self.parameters['mass']
                },
                'diameter': {
                    '_default': self.parameters['diameter'],
                },
                'velocity': {
                    '_default': self.parameters['PD1n_migration'],
                    '_updater': 'set',
                    '_divider': {
                        'divider': set_velocity_default,
                         },
                    '_emit': True,
                },
                'exchange': {
                    'IFNg': {
                        '_default': 0,
                        '_updater': 'accumulate',
                        '_divider': 'split',
                    }},
                'external': {
                    mol_id: {
                        '_default': 0.0 * CONCENTRATION_UNIT,
                        '_emit': True,
                    } for mol_id in self.parameters['external_molecules']
                },
            },
            'neighbors': {
                'present': {
                    'PD1': {
                        '_default': 0,
                        '_updater': 'set',
                    },  # membrane protein, promotes T-cell death
                    'TCR': {
                        '_default': 0,
                        '_emit': True,
                        '_updater': 'set',
                    }  # level of TCR that interacts with MHCI on tumor
                },
                'accept': {
                    'PDL1': {
                        '_default': 0,
                    },
                    'MHCI': {
                        '_default': 0,
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
        TCR_timer = states['internal']['TCR_timer']
        velocity_timer = states['internal']['velocity_timer']
        TCR = states['neighbors']['present']['TCR']
        refractory_count = states['internal']['refractory_count']
        PD1n_divide_counts = states['globals']['PD1n_divide_count']

        # death
        if cell_state == 'PD1n':
            prob_death = get_probability_timestep(
                self.parameters['death_PD1n_14hr'],
                50400,  # 14 hours (14*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_death:
                return {
                    'globals': {
                        'death': 'PD1n_apoptosis'}}

        elif cell_state == 'PD1p':
            if PDL1 >= self.parameters['PDL1_critical_number']:
                prob_death = get_probability_timestep(
                    self.parameters['death_PD1p_next_to_PDL1p_14hr'],
                    50400,  # 14 hours (14*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_death:
                    return {
                        'globals': {
                            'death': 'PD1p_PDL1_death'}}

            else:
                prob_death = get_probability_timestep(
                    self.parameters['death_PD1p_14hr'],
                    50400,  # 14 hours (14*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_death:
                    return {
                        'globals': {
                            'death': 'PD1p_apoptosis'}}

        # division
        elif cell_state == 'PD1n':
            prob_divide = get_probability_timestep(
                self.parameters['PD1n_growth_28hr'],
                100800,  # 28 hours (28*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_divide:
                return {
                    'globals': {
                        'divide': True,
                    }
                }

        elif cell_state == 'PD1p':
            prob_divide = get_probability_timestep(
                self.parameters['PD1p_growth_28hr'],
                100800,  # 28 hours (28*60*60 seconds)
                timestep)
            if random.uniform(0, 1) < prob_divide:
                PD1p_divide_count = 1
                return {
                    'globals': {
                        'divide': True,
                        'PD1p_divide_count': PD1p_divide_count}}

        elif cell_state == 'delay':
            prob_divide = probability_of_occurrence_within_interval(
                timestep, self.parameters['LymphNode_delay_growth'])
            if random.uniform(0, 1) < prob_divide:
                return {
                    'globals': {
                        'divide': True,
                    }
                }

        ## Build up an update
        update = {
            'internal': {},
            'boundary': {},
            'neighbors': {'present': {}, 'accept': {}, 'transfer': {}}}

        # TCR downregulation after 6 hours of activation
        if TCR_timer > self.parameters['activation_time']:
            TCR = self.parameters['TCR_downregulated']
            update['neighbors']['present'].update({
                'TCR': TCR})

        # update the TCR_timer for each timestep in contact with MHCI or during refractory period regardless
        if MHCI > self.parameters['ligand_threshold']:
            update['internal'].update({
                'TCR_timer': timestep})
        elif TCR_timer > self.parameters['activation_time']:
            update['internal'].update({
                'TCR_timer': timestep})

        # After the refractory period is over then
        if TCR_timer > self.parameters['activation_refractory_time']:
            TCR = self.parameters['TCR_upregulated']
            update['internal'].update({
                'TCR_timer': -self.parameters['activation_refractory_time'],
                'total_cytotoxic_packets': -states['internal']['total_cytotoxic_packets']})
            update['neighbors']['present'].update({
                'TCR': TCR})
            refractory_count_add = 1
            update['internal'].update({
                'refractory_count': refractory_count_add})

        # state transition
        new_cell_state = cell_state
        if cell_state == 'PD1n':
            if refractory_count > self.parameters['refractory_count_threshold'] or\
                    PD1n_divide_counts > self.parameters['PD1n_divide_threshold']:
                #print('PD1n become PD1p!')
                new_cell_state = 'PD1p'
                cell_state_count = 1
                update['internal'].update({
                    'cell_state': new_cell_state,
                    'cell_state_count': cell_state_count,})

        elif cell_state == 'PD1p':
            cell_state_count = 0
            update['internal'].update({
                'cell_state': new_cell_state,
                'cell_state_count': cell_state_count})

        # behavior
        # Killing by passing cytotoxic packets to tumor
        if new_cell_state == 'PD1n':

            # Update velocity timer
            if MHCI > 0:
                update['internal'].update({
                    'velocity_timer': timestep})

            # check conditional 4 fold reduction
            if MHCI >= self.parameters['ligand_threshold'] and TCR >= self.parameters['ligand_threshold']:

                # Max production for both happens for PD1- T cells in contact with MHCI+ tumor
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1n_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg
                IFNg = self.parameters['PD1n_IFNg_production'] * timestep
                update['boundary'].update({
                    'exchange': {'IFNg': int(IFNg)}})

            elif MHCI > 0 and TCR >= self.parameters['ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1n_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] / self.parameters[
                        'MHCIn_reduction_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg
                IFNg = self.parameters['PD1n_IFNg_production'] / \
                       self.parameters['MHCIn_reduction_production'] * timestep
                update['boundary'].update({
                    'exchange': {'IFNg': int(IFNg)}})

            # TODO (ERAN) -- safe to remove?
            # if states['globals']['LN_no_migration']:
            #     update['boundary'].update({
            #         'velocity': self.parameters['migration_MHCIp_tumor_dwell_velocity']})

            # Reset the velocity timer after refractory period
            elif velocity_timer >= self.parameters['PD1n_migration_refractory_time']:
                update['internal'].update({
                    'velocity_timer': {
                        '_updater': 'set',
                        '_value': 0
                    }
                })
                update['boundary'].update({
                    'velocity': self.parameters['PD1n_migration']
                })
            elif velocity_timer > self.parameters['PD1n_migration_MHCIp_tumor_dwell_time']:
                update['boundary'].update({
                    'velocity': self.parameters['PD1n_migration']})
            elif velocity_timer > 0 and velocity_timer < self.parameters['PD1n_migration_MHCIp_tumor_dwell_time']:
                update['boundary'].update({
                    'velocity': self.parameters['migration_MHCIp_tumor_dwell_velocity']})
            elif velocity_timer == 0:
                update['boundary'].update({
                    'velocity': self.parameters['PD1n_migration']})

        elif new_cell_state == 'PD1p':
            PD1 = self.parameters['PD1p_PD1_equilibrium']
            update['neighbors']['present'].update({
                'PD1': PD1})

            # Update velocity timer
            if MHCI > 0:
                update['internal'].update({
                    'velocity_timer': timestep})

            if MHCI >= self.parameters['ligand_threshold'] and TCR >= self.parameters['ligand_threshold']:

                # Max production for both happens for T cells in contact with MHCI+ tumor
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1p_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg
                IFNg = self.parameters['PD1p_IFNg_production'] * timestep
                update['boundary'].update({
                    'exchange': {'IFNg': int(IFNg)}})

            elif MHCI > 0 and TCR >= self.parameters['ligand_threshold']:

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                if states['internal']['total_cytotoxic_packets'] < self.parameters['PD1p_cytotoxic_packets_max']:
                    new_cytotoxic_packets = self.parameters['cytotoxic_packet_production'] / \
                                        self.parameters['MHCIn_reduction_production'] * timestep
                    update['internal'].update({
                        'total_cytotoxic_packets': new_cytotoxic_packets})

                # produce IFNg
                IFNg = self.parameters['PD1p_IFNg_production'] / \
                       self.parameters['MHCIn_reduction_production'] * timestep
                update['boundary'].update({
                    'exchange': {'IFNg': int(IFNg)}})

            # if states['globals']['LN_no_migration']:
            #     update['boundary'].update({
            #         'velocity': self.parameters['migration_MHCIp_tumor_dwell_velocity']})

            # Reset the velocity timer after refractory period
            elif velocity_timer >= self.parameters['PD1p_migration_refractory_time']:
                update['internal'].update({
                    'velocity_timer': {
                        '_updater': 'set',
                        '_value': 0
                    }
                })
                update['boundary'].update({
                    'velocity': self.parameters['PD1p_migration']
                })
            elif velocity_timer > self.parameters['PD1p_migration_MHCIp_tumor_dwell_time']:
                update['boundary'].update({
                    'velocity': self.parameters['PD1p_migration']})
            elif velocity_timer > 0 and velocity_timer < self.parameters['PD1p_migration_MHCIp_tumor_dwell_time']:
                update['boundary'].update({
                    'velocity': self.parameters['migration_MHCIp_tumor_dwell_velocity']})
            elif velocity_timer == 0:
                update['boundary'].update({
                    'velocity': self.parameters['PD1p_migration']})

        # keep cytotoxic transfer to max limit between tumor and T cells and remove from T cell total when transfer
        max_cytotoxic_transfer = self.parameters['cytotoxic_transfer_rate'] - states['neighbors']['transfer']['cytotoxic_packets']
        cytotoxic_transfer = min(max_cytotoxic_transfer,states['internal']['total_cytotoxic_packets'])
        update['internal'].update({
            'cytotoxic_packets': -cytotoxic_transfer})
        update['neighbors']['transfer'].update({
            'cytotoxic_packets': cytotoxic_transfer})

        return update


# test functions
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
    timeseries = simulate_process(t_cell_process, settings)

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
    out_dir='out'
):

    tcell_override = {
        '_schema': {
            'globals': {
                'death': {'_emit': False}},
            'internal': {
                'cell_state_count': {'_emit': False},
                'cell_state': {'_emit': False},
                'refractory_count': {'_emit': True},
                'total_cytotoxic_packets': {'_emit': True},
                'TCR_timer': {'_emit': True},
                'velocity_timer': {'_emit': False}},
            'neighbors': {
                'present': {
                    'PD1': {'_emit': True}},
                'accept': {
                    'PDL1': {'_emit': True},
                    'MHCI': {'_emit': True}}},
            'boundary': {
                'velocity': {'_emit': False},
                'external': {
                    'IFNg': {'_emit': False}}}}}

    combined_raw_data = {}
    for single_idx in range(batch_size):
        t_cell_process = TCellProcess(tcell_override)
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
        raw_data = simulate_process(t_cell_process, sim_settings)
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

    parser = argparse.ArgumentParser(description='run t_cell process')
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
            timeline=timeline,
            out_dir=out_dir)
