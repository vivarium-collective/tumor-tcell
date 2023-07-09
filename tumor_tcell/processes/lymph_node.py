"""
Lymph Node Environment Process
"""

import math
import random
from vivarium.core.process import Process
from vivarium.core.engine import pp
from tumor_tcell.processes.t_cell import get_probability_timestep, TIMESTEP
from vivarium.library.units import remove_units
from tumor_tcell.library.location import random_location, DEFAULT_BOUNDS


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
        'tumor_env_bounds': remove_units(DEFAULT_BOUNDS),
        'time_dendritic_finds_tcell': 2.0,  # hours TODO - I don't think we use this
        'n_tcells_in_lymph_node': 3,  # number of cells in 3D LN structure. CD8+ T cells 12% of all cells, \
            # 5x10^6 cells/LN, so 600,000 T cells. 0.5% of those are antigen specific, so about 3000 T cells \
            # within the lymph node that are reactive (total pool of cells that could proliferated/divide. \
            # - from data and divide 3000/1000 (1/1000th of space that we are simulating) = 3
        'tcell_find_dendritic_time': 0.95,  # 95% will find dendritic in 4 hrs (Itano, 2003);;(Bousso, 2008)
        'expected_dendritic_transit_time': 20,  # 28800 8*60*60. 8 hour delay between the time that a dendritic \ TODO - @John change back
            # cell leaves microenvironment until it is ready to interact with t cells in the LN and interact with \
            # T cells that take about 4 hours to find it for a total of 12 hours total until engagement is \
            # seen (Itano, 2003);;(Bousso, 2008)
        'expected_tcell_transit_time': 30,  # 3600 60*60. arrive in tumor environment after 1 hour of migration in,\ TODO - @John change back
            # efferent lymph to circulation (Hunter, 2016)
        'expected_division_interval': 10,  # 14400 divide approximately every 4 hours, or 5-6 times in 24 hours. \ TODO - @John change back
            # 3*60*60=10800, (Mempel, 2004);(Bousso, 2008)
        'expected_interaction_duration': 20,  # 28800 8*60*60 t cells interact with dendritic cells for approximately \ TODO - @John change back
            # 8 hours (Itano, 2003)
        'expected_delay_before_migration': 40,  # 43200 12*60*60. t cells wait approx 12 hours after interaction is \ TODO - @John change back
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
        tumor_env_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'}},
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {},
                    'location': {}
                },
        }}

        lymph_node_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'},
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {
                        '_default': '',
                    }}}}

        # TODO -- reuse the schemas more instead of copying
        in_transit_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'},
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {
                        '_default': '',
                    }},
        }}

        return {
            'cells': tumor_env_schema,
            'lymph_node': lymph_node_schema,
            'in_transit': in_transit_schema
        }

    def next_update(self, timestep, states):
        microenvironment_cells = states['cells']
        lymph_node_cells = states['lymph_node']
        in_transit = states['in_transit']

        update = {
            'cells': {
                '_add': [],
                '_delete': [],
            },
            'lymph_node': {
                '_add': [],
                '_delete': [],
            },
            'in_transit': {
                '_add': [],
                '_delete': [],
            }
        }

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
                    # interact with dendritic cells for 8 hours, then 12 hour activation delay, and then migrate to tumor
                    prob_interaction_completion = probability_of_occurrence_within_interval(
                        timestep, self.parameters['expected_interaction_duration'])
                    if random.uniform(0, 1) < prob_interaction_completion:
                        # first delay, then migrate
                        update['lymph_node'][cell_id] = {
                            'internal': {'cell_state': 'delay'}
                        }

                elif cell_state == 'delay':
                    # get probability of migration starting
                    prob_migration = probability_of_occurrence_within_interval(
                        timestep, self.parameters['expected_delay_before_migration'])
                    if random.uniform(0, 1) < prob_migration:
                        # TODO -- move it to "in transit", which should be relatively fast
                        # begin transit
                        update['in_transit']['_add'].append({'key': cell_id, 'state': specs})
                        update['cells']['_delete'].append(cell_id)
                    else:
                        # start dividing over next 12-18 hours. Probability of division goes up
                        prob_divide = probability_of_occurrence_within_interval(
                            timestep, self.parameters['expected_division_interval'])
                        if random.uniform(0, 1) < prob_divide:
                            update['lymph_node'][cell_id].update({
                                'globals': {
                                    'divide': True}})
                else:
                    # Calculate probability of finding/initializing interaction with dendritic cells
                    # TODO -- this should depend on dendritic cell being present. Not interacting alone
                    prob_interaction = get_probability_timestep(
                        self.parameters['tcell_find_dendritic_time'],
                        10,  # 14400 6 hours (6*60*60 seconds) TODO - @John change back
                        timestep) #(Itano, 2003)
                    if random.uniform(0, 1) < prob_interaction:
                        # this t-cell is now interacting
                        update['lymph_node'][cell_id] = {
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
                    # begin transit
                    update['in_transit']['_add'].append({'key': cell_id, 'state': specs})
                    update['cells']['_delete'].append(cell_id)

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
                    # move to lymph node
                    update['lymph_node']['_add'].append({'key': cell_id, 'state': specs})
                    update['in_transit']['_delete'].append(cell_id)
            if cell_type == 't-cell':
                # t cells move from LN to tumor
                prob_arrival = probability_of_occurrence_within_interval(
                    timestep, self.parameters['expected_tcell_transit_time'])
                if random.uniform(0, 1) < prob_arrival:
                    # move to lymph node
                    location = random_location(self.parameters['tumor_env_bounds'])
                    specs['boundary']['location'] = location
                    update['cells']['_add'].append({'key': cell_id, 'state': specs})
                    update['in_transit']['_delete'].append(cell_id)


        # print(f'STATES TUMOR ENV: {list(states["cells"].keys())}')
        # print(f'STATES LN: {list(states["lymph_node"].keys())}')
        # print(f'STATES IN TRANSIT: {list(states["in_transit"].keys())}')
        # print(f'UPDATE TUMOR ENV: {update["cells"]}')
        # print(f'UPDATE LN: {update["lymph_node"]}')
        # print(f'UPDATE IN TRANSIT: {update["in_transit"]}')
        if update['in_transit'].get('_add') or update['in_transit'].get('_delete'):
            print(f'UPDATE: {pp(update)}')
            x=0

        return update


from vivarium.core.engine import Engine
def test_lymph_node():
    n_in_cells = 4
    n_in_transit = 2
    n_in_ln = 1
    simtime = 1000000

    cell_schema = {
        'internal': {
            'cell_state': 'inactive'},
        'boundary': {
            'cell_type': '',
            'location': {}}
    }

    ln = LymphNode()

    print(pp(ln.ports_schema()))

    sim = Engine(
        processes={'ln': ln},
        topology={'ln': {
            'cells': ('cells',),
            'lymph_node': ('lymph_node',),
            'in_transit': ('in_transit',)}
        },
        initial_state={
            'cells': {f'c{i}': cell_schema for i in range(n_in_cells)},
            'lymph_node': {f'ln{i}': cell_schema for i in range(n_in_ln)},
            'in_transit': {f't{i}': cell_schema for i in range(n_in_transit)},
        }
    )

    sim.update(simtime)

    print(pp(sim.state.get_value()))
    x=0



if __name__ == '__main__':
    test_lymph_node()
