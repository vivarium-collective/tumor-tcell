"""
======================
Tumor microenvironment
======================
"""

from vivarium.core.process import Generator
from tumor_tcell.processes.multibody_neighbors import MultibodyNeighbors



class TumorMicroEnvironment(Generator):
    """ Tumor microenvironment compartment
    Models the environment in which tumors grow
    """
    defaults = {
        'multibody': {
            'bounds': [10, 10]  # the size of the microenvironment
        },
        'chemical_fields': {},
        '_schema': {},
    }

    def __init__(self, config=None):
        super(TumorMicroEnvironment, self).__init__(config)

    def generate_processes(self, config):
        multibody = MultibodyNeighbors(config['multibody'])
        return {
            'multibody': multibody
        }

    def generate_topology(self, config):
        return {
            'multibody': {
                'cells': ('cells',)
            }
        }

