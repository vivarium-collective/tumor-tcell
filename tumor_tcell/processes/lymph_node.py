"""
Lymph Node Environment Process
"""


from vivarium.core.process import Process
from tumor_tcell.processes.neighbors import daughter_locations


class LymphNode(Process):
    """Enable exchange of Dendritic cells between location,
    interaction of T cells and Dendritic cells,
    and T cells leaving to the environment
    """
    defaults = {
        'time_dendritic_finds_tcell': 2.0,  # hours
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        cell_schema = {
            '*': {
                'internal': {
                    'cell_state': {
                        '_default': 'inactive',
                        '_emit': True,
                        '_updater': 'set',
                    },
                },
                'boundary': {
                    # cell_type must be either 'tumor', 't_cell', or 'dendritic'
                    'cell_type': {},
                    # 'location': {
                    #     '_emit': True,
                    #     '_default': [
                    #         0.5 * bound for bound in self.parameters['bounds']],
                    #     '_updater': 'set',
                    #     '_divider': {
                    #         'divider': daughter_locations,
                    #         'topology': {
                    #             'diameter': ('..', 'diameter',)
                    #         }
                    #     }
                    # },
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
                        '*': {'_default': 0.0}}}}}

        return {
            'cells': cell_schema,
            'lymph_node': {},
        }

    def next_update(self, timestep, states):
        microenvironment_cells = states['cells']
        lymph_node_cells = states['lymph_node']

        # get the cells that are in the tumor environment

        # get the cells that are in the lymph node

        update = {}
        # move cells between tumor environment and lymph node


        return update


def test_lymph_node():
    pass


if __name__ == '__main__':
    test_lymph_node()
