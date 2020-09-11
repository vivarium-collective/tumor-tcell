
from vivarium.core.process import Generator
from vivarium.processes.meta_division import MetaDivision
from tumor_tcell.processes.tumor import TumorProcess


class TumorCompartment(Generator):

    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        '_schema': {},
    }

    def __init__(self, config):
        super(TumorCompartment, self).__init__(config)

    def generate_processes(self, config):
        daughter_path = config['daughter_path']
        agent_id = config['agent_id']

        division_config = dict(
            config.get('division', {}),
            daughter_path=daughter_path,
            agent_id=agent_id,
            compartment=self)

        Tumor = TumorProcess(config.get('growth', {}))
        division = MetaDivision(division_config)

        return {
            'Tumor': Tumor,
            'division': division}

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        agents_path = config['agents_path']
        return {
            'Tumor': {
                'internal': ('internal',),
                'boundary': boundary_path,
                'global': boundary_path},
            'division': {
                'global': boundary_path,
                'cells': agents_path},
            }