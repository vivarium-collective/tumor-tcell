"""
======================
Tumor microenvironment
======================
"""

from vivarium.core.process import Generator
from tumor_tcell.processes.neighbors import Neighbors
from tumor_tcell.processes.fields import Fields


class TumorMicroEnvironment(Generator):
    """ Tumor microenvironment compartment
    Models the environment in which tumors grow
    """
    defaults = {
        'neighbors': {
            'bounds': [10, 10]  # the size of the microenvironment
        },
        'fields': {},
        '_schema': {},
    }

    def __init__(self, config=None):
        super(TumorMicroEnvironment, self).__init__(config)

    def generate_processes(self, config):

        # initialize processes
        neighbors = Neighbors(config['neighbors'])

        # make dictionary of processes
        return {
            'neighbors': neighbors
        }

    def generate_topology(self, config):
        return {
            'neighbors': {
                'cells': ('cells',)
            }
        }
