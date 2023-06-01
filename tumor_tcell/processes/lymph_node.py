"""
Lymph Node Environment Process
"""


from vivarium.core.process import Process
from vivarium.core.engine import pp


class LymphNode(Process):
    """Enable exchange of Dendritic cells between location,
    interaction of T cells and Dendritic cells,
    and T cells leaving to the environment
    """
    defaults = {
        'time_dendritic_finds_tcell': 2.0,  # hours
        'n_tcells_in_lymph_node': 100,  # number of cells in 3D LN structure. TODO -- calculate number in a 2D slice
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def initial_state(self, config):
        agent_ids = config.get('agent_ids')
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
                        '_updater': 'set',
                    },
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {},
                    'location': {},  # TODO -- needed?
                },
                'neighbors': {
                    'present': {
                        '*': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'accept': {
                        '*': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'transfer': {
                        '*': {'_default': 0.0}},
                    'receive': {
                        '*': {'_default': 0.0}
                    }
                }}}

        lymph_node_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set',
                    },
                    'lymph_node_timer': {
                        '_default': 0
                    }  # TODO -- maybe this should just be for T cells rather than everything?
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {},
                    'location': {},  # TODO -- needed?
                },
            }
        }

        return {
            'cells': cell_schema,
            'lymph_node': lymph_node_schema,
        }

    def next_update(self, timestep, states):
        microenvironment_cells = states['cells']
        lymph_node_cells = states['lymph_node']

        update = {
            'cells': {
                '_add': [],
                '_remove': [],
            },
            'lymph_node': {
                '_add': [],
                '_remove': [],
            }
        }

        #############################
        # interactions and division #
        #############################

        # interactions between dendritic and t cells
        # This will be about setting timers on the t cells and tracking interaction
        # TODO if dendritic cells are in lymph node, then start random timer, to set time for T cells to get the dendritic cell -- different threshold? different timer?
        # if interaction, increase the T cell lymph_node_interaction_timer state

        # TODO -- do this in the T cell process: let them divide in the LN
        # Maybe increase division probability 12-16 hours of division.


        ################
        # moving cells #
        ################

        # move dendritic cells from tumor to the lymph node
        # Once dendritic cells are active, they move to LN and stay there until chemokines subsist. We are assuming they remain in LN for the duration of the simulation
        for cell_id, specs in microenvironment_cells.items():
            if specs['boundary']['cell_type'] == 'dendritic':
                # get the state of the dendritic cells in the microenvironment
                if specs['internal']['cell_state'] == 'active': # either active or inactive
                    update['lymph_node']['_add'].append({cell_id: specs})  # if active, move to lymph node, remove from microenvironment
                    update['cells']['_remove'].append({cell_id: specs})
                    # 12 hour delay between the time that a dendritic cell leaves microenvironment until it is ready to interact with t cells in the LN

        # move T cells from lymph node to the tumor
        # T cells move one way from the LN to the tumor. going back to LN is more rare, so we are leaving it out
        for cell_id, specs in lymph_node_cells.item():
            if specs['boundary']['cell_type'] == 't-cell':
                pass
                # use timers for the t cells.
                # We don't want t cells to become refractory in the LN.

                # T cells interact with dendritic cells for 8 hours (on average -- get a distribution), then 12 hour activation delay, another hour of interactions, and then migrate to tumor
                # Random timer for each T cell to



        return update


def test_lymph_node():
    pass


if __name__ == '__main__':
    test_lymph_node()
