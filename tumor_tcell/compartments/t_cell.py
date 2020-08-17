
from vivarium.core.process import Generator
from tumor_tcell.processes.t_cell import TCellProcess


class TCellCompartment(Generator):

    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple()}

    def __init__(self, config):
        self.config = config
        for key, value in self.defaults.items():
            if key not in self.config:
                self.config[key] = value

        # paths
        self.boundary_path = config.get('boundary_path', self.defaults['boundary_path'])
        self.agents_path = config.get('agents_path', self.defaults['agents_path'])

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
        return {
            't_cell': {
                'internal': ('internal',),
                'boundary': self.boundary_path,
                'global': self.boundary_path},
            'division': {
                'global': self.boundary_path,
                'cells': self.agents_path},
            }