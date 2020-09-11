
from vivarium.core.process import Generator
from vivarium.processes.meta_division import MetaDivision
from tumor_tcell.processes.t_cell import TCellProcess


class TCellCompartment(Generator):

    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        '_schema': {},
    }

    def __init__(self, config):
        super(TCellCompartment, self).__init__(config)

    def generate_processes(self, config):
        daughter_path = config['daughter_path']
        agent_id = config['agent_id']

        division_config = dict(
            config.get('division', {}),
            daughter_path=daughter_path,
            agent_id=agent_id,
            compartment=self)

        t_cell = TCellProcess(config.get('growth', {}))
        division = MetaDivision(division_config)

        return {
            't_cell': t_cell,
            'division': division}

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']
        return {
            't_cell': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'global': boundary_path},
            'division': {
                'global': boundary_path,
                'cells': agents_path},
            }