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
                        '_updater': 'set' }},
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {}},
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
                        '*': {'_default': 0.0}}}}}

        lymph_node_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_updater': 'set'},
                    'lymph_node_timer': {
                        '_default': 0}},
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {}}}}

        return {
            'cells': cell_schema,
            'lymph_node': lymph_node_schema}

    def next_update(self, timestep, states):
        microenvironment_cells = states['cells']
        lymph_node_cells = states['lymph_node']

        update = {
            'cells': {
                '_add': [],
                '_delete': [],
            },
            'lymph_node': {
                '_add': [],
                '_delete': [],
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


        ##############
        # lymph node #
        ##############

        # move T cells from lymph node to the tumor
        # T cells move one way from the LN to the tumor. going back to LN is more rare, so we are leaving it out

        for cell_id, specs in lymph_node_cells.items():
            cell_type = specs['boundary']['cell_type']
            cell_state = specs['internal']['cell_state']

            if cell_type == 't-cell':
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
                    # if active, move to lymph node, remove from microenvironment
                    update['lymph_node']['_add'].append({'key': cell_id, 'state': specs})
                    update['cells']['_delete'].append(cell_id)
                elif cell_state == 'inactive':
                    update['cells'][cell_id] = {
                        'internal': {
                            'lymph_node_timer': timestep  # increase the timer
                        }
                    }

                    # 12 hour delay between the time that a dendritic cell leaves microenvironment until it is ready to interact with t cells in the LN


        print(f'STATES CELLS: {list(states["cells"].keys())}')
        print(f'STATES LN: {list(states["lymph_node"].keys())}')
        print(f'UPDATE CELLS: {update["cells"]}')
        print(f'UPDATE LN: {update["lymph_node"]}')

        return update


def test_lymph_node():
    pass


if __name__ == '__main__':
    test_lymph_node()
