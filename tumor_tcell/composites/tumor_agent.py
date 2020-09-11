
from vivarium.core.process import Generator
from vivarium.processes.meta_division import MetaDivision
from tumor_tcell.processes.tumor import TumorProcess


class TumorAgent(Generator):

    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        'divide': True,
        '_schema': {},
    }

    def __init__(self, config):
        super(TumorAgent, self).__init__(config)

    def generate_processes(self, config):

        # initialize processes
        tumor = TumorProcess(config.get('tumor', {}))

        # make dictionary of processes
        processes = {
            'tumor': tumor,
        }

        # if divide set to true, add meta-division processes
        if config['divide']:
            daughter_path = config['daughter_path']
            agent_id = config['agent_id']
            meta_division_config = dict(
                {},
                daughter_path=daughter_path,
                agent_id=agent_id,
                compartment=self)
            meta_division = MetaDivision(meta_division_config)
            processes['division'] = meta_division

        return processes


    def generate_topology(self, config):

        # retrieve paths
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']

        # make topology by mapping ports
        topology = {
            'tumor': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'globals': boundary_path,
                'neighbors': ('neighbors',),
            },
        }
        if config['divide']:
            topology.update({
                'division': {
                    'global': boundary_path,
                    'agents': agents_path,
                }})
        return topology
