"""
Lymph Node Environment Process
"""

import random
from vivarium.core.process import Process
from vivarium.core.engine import pp
from tumor_tcell.processes.t_cell import get_probability_timestep


class LymphNode(Process):
    """Enable exchange of Dendritic cells between location,
    interaction of T cells and Dendritic cells,
    and T cells leaving to the environment
    """
    defaults = {
        'time_dendritic_finds_tcell': 2.0,  # hours
        'n_tcells_in_lymph_node': 100,  # number of cells in 3D LN structure. TODO -- calculate number in a 2D slice
        'transit_time': 12.0,  # 12 hour delay between the time that a dendritic cell leaves microenvironment until it is ready to interact with t cells in the LN
        'tcell_find_dendritic_6hr': 0.5,  # TODO(ERAN): 0.95 will find dendritic in 6 hrs (actually 6-8)
        'prob_division_when_stimulated': 0.20,  # PD1p_growth_28hr = 20% division in 28 hours. TODO (ERAN) recalculate. Should divide 4-6 times in 12-16 hours
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def initial_state(self, config=None):
        if config:
            agent_ids = config.get('agent_ids', [])
        # TODO all of these need to be of 'cell_state': 'PD1n'
        # TODO -- T cells need a random timer value between 0 and 8 hours to start with
        return {
            'lymph_node': {}  # TODO put tcells in lymph node (maybe in main.py)
        }

    def ports_schema(self):
        cell_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'}},
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {}},
        }}

        lymph_node_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'},
                    'lymph_node_timer': {
                        '_default': 0
                    }
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {}}}}

        # TODO -- reuse the schemas more instead of copying
        in_transit_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'},
                    'transit_timer': {
                        '_default': 0
                    }
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {}},
        }}

        return {
            'cells': cell_schema,
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

        # interactions between dendritic and t cells
        # This will be about setting timers on the t cells and tracking interaction
        # TODO if dendritic cells are in lymph node, then start random timer, to set time for T cells to get the dendritic cell -- different threshold? different timer?
        # if interaction, increase the T cell lymph_node_interaction_timer state

        # TODO -- do this in the T cell process: let them divide in the LN
        # Maybe increase division probability 12-16 hours of division.


        ##############
        # lymph node #
        ##############

        # move T cells from lymph node to the tumor
        # T cells move one way from the LN to the tumor. going back to LN is more rare, so we are leaving it out
        # ~60k-180k total T cells. 0.005 are antigen-specific. We want a 2D slice ~1
        # Start off with ~3 t cells in LN, allow them to interact

        for cell_id, specs in lymph_node_cells.items():
            cell_type = specs['boundary']['cell_type']
            cell_state = specs['internal']['cell_state']

            if cell_type == 'dendritic':
                pass
                # update['lymph_node'][cell_id] = {
                #     'internal': {
                #         'lymph_node_timer': timestep  # increase the timer
                #     }
                # }
            if cell_type == 't-cell':
                # Sense dendritic cells. We are assuming all t cells are a match to the dendritic cells.
                # Calculate probability of interaction
                prob_interaction = get_probability_timestep(
                    self.parameters['tcell_find_dendritic_6hr'],
                    21600,  # 6 hours (6*60*60 seconds)
                    timestep)
                if random.uniform(0, 1) < prob_interaction:
                    update['lymph_node'][cell_id] = {}
                    # this t-cell is now interacting
                    # start dividing over next 12-18 hours. Probability of division goes up
                    prob_divide = get_probability_timestep(
                        self.parameters['prob_division_when_stimulated'],
                        100800,  # 28 hours (28*60*60 seconds) TODO ???
                        timestep)
                    if random.uniform(0, 1) < prob_divide:
                        update['lymph_node'][cell_id].update({'globals': {'divide': True}})

                    # time how long since interaction
                    update['lymph_node'][cell_id].update({
                        'internal': {
                            'lymph_node_timer': timestep  # increase the timer
                        }
                    })


                # in LN, t-cells are confined to different spaces "t cell zones"
                # during immune response, dendritic cells move to these regions.  T cells sense dendritic cells, and if specific match there are repeated interactions for ~6-8 hours. These cells are acitvated and start dividing (4 divisions in next 2 days). and then leav to tumor.

                if cell_state == 'PD1p':
                    pass

                # use timers for the t cells.
                # We don't want t cells to become refractory in the LN.

                # T cells interact with dendritic cells for 8 hours (on average -- get a distribution), then 12 hour activation delay, another hour of interactions, and then migrate to tumor
                # Random timer for each T cell to


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
                if specs['internal']['transit_timer'] >= self.parameter['transit_time']:  # TODO dendritic_transit_timer
                    specs['internal']['transit_timer'] = 0.0
                    # move to lymph node
                    update['lymph_node']['_add'].append({'key': cell_id, 'state': specs})
                    update['in_transit']['_delete'].append(cell_id)
                else:
                    update['lymph_node'][cell_id] = {
                        'internal': {
                            'transit_timer': timestep  # increase the timer
                        }
                    }
            if cell_type == 't-cell':
                # t cells move from LN to tumor
                pass


        print(f'STATES TUMOR ENV: {list(states["cells"].keys())}')
        print(f'STATES LN: {list(states["lymph_node"].keys())}')
        print(f'STATES IN TRANSIT: {list(states["in_transit"].keys())}')
        print(f'UPDATE TUMOR ENV: {update["cells"]}')
        print(f'UPDATE LN: {update["lymph_node"]}')
        print(f'UPDATE IN TRANSIT: {update["in_transit"]}')

        return update


def test_lymph_node():
    pass


if __name__ == '__main__':
    test_lymph_node()
