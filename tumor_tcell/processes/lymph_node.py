"""
Lymph Node Environment Process
"""

import math
import random
from vivarium.core.process import Process
from vivarium.core.engine import pp, Engine
from tumor_tcell.processes.t_cell import get_probability_timestep, TIMESTEP
from vivarium.library.units import remove_units
from tumor_tcell.library.location import random_location, DEFAULT_BOUNDS
from tumor_tcell.processes.local_field import LENGTH_UNIT
from tumor_tcell.processes.neighbors import DEFAULT_MASS_UNIT, DEFAULT_VELOCITY_UNIT

DEFAULT_CELL_TYPE = 'default_cell_type'


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


class LymphNode(Process):
    """Enable exchange of Dendritic cells between location,
    interaction of T cells and Dendritic cells,
    and T cells leaving to the environment
    """
    defaults = {
        'time_step': TIMESTEP,
        'tumor_env_bounds': DEFAULT_BOUNDS,
        'time_dendritic_finds_tcell': 2.0,  # hours TODO - I don't think we use this
        'n_tcells_in_lymph_node': 3,  # number of cells in 3D LN structure. CD8+ T cells 12% of all cells, \
            # 5x10^6 cells/LN, so 600,000 T cells. 0.5% of those are antigen specific, so about 3000 T cells \
            # within the lymph node that are reactive (total pool of cells that could proliferated/divide. \
            # - from data and divide 3000/1000 (1/1000th of space that we are simulating) = 3
        'tcell_find_dendritic_time': 0.95,  # 95% will find dendritic in 4 hrs (Itano, 2003);;(Bousso, 2008)
        'expected_dendritic_transit_time': 60,  # Assuming that some DCs already present within lymph node 28800 8*60*60. \
            # 8 hour delay between the time that a dendritic \
            # cell leaves microenvironment until it is ready to interact with t cells in the LN and interact with \
            # T cells that take about 4 hours to find it for a total of 12 hours total until engagement is \
            # seen (Itano, 2003);;(Bousso, 2008)
        'expected_tcell_transit_time': 3600,  # 3600 60*60. arrive in tumor environment after 1 hour of migration in,\
            # efferent lymph to circulation (Hunter, 2016)
        # 'expected_division_interval': 14400,  # 14400 divide approximately every 4 hours, or 5-6 times in 24 hours. \
        #     # 3*60*60=10800, (Mempel, 2004);(Bousso, 2008)
        'expected_interaction_duration': 28800,  # 28800 8*60*60 t cells interact with dendritic cells for approximately \
            # 8 hours (Itano, 2003)
        'expected_delay_before_migration': 43200,  # 43200 12*60*60. t cells wait approx 12 hours after interaction is \
            # complete before starting migration (Itano, 2003);(Bousso, 2008)
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def initial_state(self, config=None):
        if config:
            agent_ids = config.get('agent_ids', [])
        return {
            'lymph_node': {}  # TODO put tcells in lymph node (maybe in main.py)
        }

    def ports_schema(self):

        agents_schema = {
            'agents': {
                '*': {
                    'internal': {
                        'cell_state': {
                            '_default': 'inactive',
                            '_updater': 'set',
                            '_emit': True,}},
                    'boundary': {
                        # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                        'cell_type': {'_default': DEFAULT_CELL_TYPE,
                                      '_emit': True,},
                        'diameter': {'_default': 1.0 * LENGTH_UNIT,
                                     '_emit': True,},
                        'mass': {'_default': 1.0 * DEFAULT_MASS_UNIT,
                                 '_emit': True,},
                        'velocity': {'_default': 0.0 * DEFAULT_VELOCITY_UNIT,
                                     '_emit': True,},
                        'location': {
                            '_default': [0.5 * bound for bound in self.parameters['tumor_env_bounds']],
                            '_updater': 'set',
                            '_emit': True,
                        },
                        'external': {'IFNg': {'_default': 0.0,
                                              '_emit': True,},
                                     'tumor_debris': {'_default': 0.0,
                                                      '_emit': True,},},  # TODO -- this should not be required here
                        'death': {'_default': False,
                                  '_emit': True,}  # TODO -- this should not be required here
                    },
                    # initialize the schema for neighbors so cells will have it when moving back to tumor
                    'neighbors': {
                        'present': {'*': {'_default': 0.0, '_emit': True,}},
                        'accept': {'*': {'_default': 0.0, '_emit': True,}},
                        'transfer': {'*': {'_default': 0.0, '_emit': True,}},
                        'receive': {'*': {'_default': 0.0, '_emit': True,}}
                    }
                }
            }
        }
        return {
            'cells': agents_schema,
            'lymph_node': agents_schema,
            'in_transit': agents_schema
        }

    def next_update(self, timestep, states):
        microenvironment_cells = states['cells']['agents']
        lymph_node_cells = states['lymph_node']['agents']
        in_transit = states['in_transit']['agents']

        cells_update = {}
        in_transit_update = {}
        lymph_node_update = {}

        ##############
        # lymph node #
        ##############

        # move T cells from lymph node to the tumor
        # T cells move one way from the LN to the tumor. going back to LN is more rare, so we are leaving it out
        # ~60k-180k total T cells. 0.005 are antigen-specific. We want a 2D slice ~1
        # Start off with ~3 t cells in LN, allow them to interact

        # check if there are dendritic cells present for interacting with T cells
        dendritic_cells_present = any([
            specs['boundary']['cell_type'] == 'dendritic'
            for specs in lymph_node_cells.values()])

        for cell_id, specs in lymph_node_cells.items():
            cell_type = specs['boundary']['cell_type']
            cell_state = specs['internal']['cell_state']

            if cell_type == 't-cell' and dendritic_cells_present:
                if cell_state == 'interacting':
                    # interact with dendritic cells for 8 hours, then 12 hour delay before migrating to tumor
                    prob_interaction_completion = probability_of_occurrence_within_interval(
                        timestep, self.parameters['expected_interaction_duration'])
                    if random.uniform(0, 1) < prob_interaction_completion:
                        # first delay, then migrate
                        lymph_node_update[cell_id] = {
                            'internal': {'cell_state': 'delay'}
                        }

                elif cell_state == 'delay':
                    # get probability of migration starting
                    prob_migration = probability_of_occurrence_within_interval(
                        timestep, self.parameters['expected_delay_before_migration'])
                    if random.uniform(0, 1) < prob_migration:
                        if '_move' not in lymph_node_update:
                            lymph_node_update['_move'] = []
                        # specs['internal']['cell_state'] = 'PD1n'  # TODO -- get this state passed to cell
                        # begin transit from lymph node
                        lymph_node_update['_move'].append({
                            'source': (cell_id,),
                            'target': ('in_transit', 'agents',),
                            'update': {'internal': {'cell_state': 'PD1n'}}
                        })

                else:
                    # Calculate probability of finding/initializing interaction with dendritic cells
                    # TODO -- this should depend on dendritic cell being present. Not interacting alone
                    prob_interaction = get_probability_timestep(   # TODO -- ERAN -- why not probability_of_occurrence_within_interval?
                        self.parameters['tcell_find_dendritic_time'],
                        14400,  # 14400 6 hours (6*60*60 seconds)
                        timestep) #(Itano, 2003)
                    if random.uniform(0, 1) < prob_interaction:
                        # this t-cell is now interacting
                        lymph_node_update[cell_id] = {
                            'internal': {'cell_state': 'interacting'}
                        }

        #####################
        # tumor environment #
        #####################

        # Move dendritic cells from tumor to the lymph node
        # Once dendritic cells are active, they move to LN and stay there until chemokines subsist.
        # We are assuming they remain in LN for the duration of the simulation

        for cell_id, specs in microenvironment_cells.items():
            cell_type = specs['boundary']['cell_type']
            cell_state = specs['internal']['cell_state']
            if cell_type == 'dendritic':
                if cell_state == 'active':
                    if '_move' not in cells_update:
                        cells_update['_move'] = []
                    # begin transit from tumor environment
                    cells_update['_move'].append({
                        'source': (cell_id,),
                        'target': ('in_transit', 'agents',),
                    })

        ##############
        # in transit #
        ##############

        for cell_id, specs in in_transit.items():
            cell_type = specs['boundary']['cell_type']

            if cell_type == 'dendritic':
                # dendritic cells move only from tumor to LN
                prob_arrival = probability_of_occurrence_within_interval(
                    timestep, self.parameters['expected_dendritic_transit_time'])
                if random.uniform(0, 1) < prob_arrival:
                    if '_move' not in in_transit_update:
                        in_transit_update['_move'] = []
                    # arrive at lymph node
                    in_transit_update['_move'].append({
                        'source': (cell_id,),
                        'target': ('lymph_node', 'agents',),
                    })
            if cell_type == 't-cell':
                # tcells move from in_transit to tumor
                prob_arrival = probability_of_occurrence_within_interval(
                    timestep, self.parameters['expected_tcell_transit_time'])
                if random.uniform(0, 1) < prob_arrival:
                    if '_move' not in in_transit_update:
                        in_transit_update['_move'] = []
                    # arrive at tumor
                    location = random_location(self.parameters['tumor_env_bounds'])
                    specs['boundary']['location'] = location  # TODO -- need to add this location in move
                    in_transit_update['_move'].append({
                        'source': (cell_id,),
                        'target': ('cells', 'agents'),
                        'update': {'boundary': {'location': location}}
                    })

        return {
            'cells': {'agents': cells_update},
            'lymph_node': {'agents': lymph_node_update},
            'in_transit': {'agents': in_transit_update},
        }



def test_lymph_node():
    simtime = 1000000
    init_state = {
        'cells': {
            'tcell_0': {'internal': {'cell_state': 'PD1n'}, 'boundary': {'cell_type': 't-cell', 'location': []}},
            'tcell_1': {'internal': {'cell_state': 'PD1p'}, 'boundary': {'cell_type': 't-cell', 'location': []}},
            'tumor_0': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_1': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_2': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_3': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_4': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_5': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_6': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_7': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_8': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_9': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_10': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_11': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_12': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_13': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_14': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_15': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_16': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_17': {'internal': {'cell_state': 'PDL1p'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_18': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'tumor_19': {'internal': {'cell_state': 'PDL1n'}, 'boundary': {'cell_type': 'tumor', 'location': []}},
            'dendritic_0': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_1': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_2': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_3': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_4': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_5': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_6': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_7': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_8': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_9': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_10': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_11': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_12': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_13': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_14': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_15': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_16': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_17': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_18': {'internal': {'cell_state': 'inactive'}, 'boundary': {'cell_type': 'dendritic', 'location': []}},
            'dendritic_19': {'internal': {'cell_state': 'active'}, 'boundary': {'cell_type': 'dendritic', 'location': []}}},
        'lymph_node': {
            'tcell_LN_0': {'internal': {'cell_state': 'PD1n'}, 'boundary': {'cell_type': 't-cell'}},
            'tcell_LN_1': {'internal': {'cell_state': 'PD1n'}, 'boundary': {'cell_type': 't-cell'}},
            'tcell_LN_2': {'internal': {'cell_state': 'PD1n'}, 'boundary': {'cell_type': 't-cell'}}},
        'in_transit': {}}

    # make the simulation and run
    ln = LymphNode()
    sim = Engine(
        processes={'ln': ln},
        topology={'ln': {
            'cells': ('cells',),
            'lymph_node': ('lymph_node',),
            'in_transit': ('in_transit',)}
        },
        initial_state=init_state
    )
    sim.update(simtime)
    data = sim.emitter.get_data()
    print(pp(data[simtime]))


if __name__ == '__main__':
    test_lymph_node()
