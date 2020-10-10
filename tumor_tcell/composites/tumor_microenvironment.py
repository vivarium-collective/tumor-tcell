"""
======================
Tumor microenvironment
======================
"""

from vivarium.core.process import Generator
from tumor_tcell.processes.neighbors import Neighbors
from vivarium_cell.processes.diffusion_field import DiffusionField


class TumorMicroEnvironment(Generator):
    """ Tumor microenvironment compartment
    Models the environment in which tumors grow
    """
    defaults = {
        'neighbors_multibody': {
            'bounds': [10, 10]  # the size of the microenvironment
        },
        'fields': {},
        '_schema': {},
    }

    def __init__(self, config=None):
        super(TumorMicroEnvironment, self).__init__(config)

    def generate_processes(self, config):

        # initialize processes
        neighbors_multibody = Neighbors(config['neighbors_multibody'])
        diffusion_fields = DiffusionField(config['fields'])

        # make dictionary of processes
        return {
            'neighbors_multibody': neighbors_multibody,
            'diffusion_fields': diffusion_fields,
        }

    def generate_topology(self, config):
        return {
            'neighbors_multibody': {
                'cells': ('agents',)
                # use agents store for integration with agent_environment_experiment in composition
                # TODO (Eran) -- update agent_environment_experiment to allow for any store name
            },
            'diffusion_fields': {
                'agents': ('agents',),
                'fields': ('fields',),
                'dimensions': ('dimensions',),
            }
        }
