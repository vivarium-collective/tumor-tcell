import copy
from vivarium.core.process import Deriver

def append_log(current_value, new_value):
    log = copy.deepcopy(current_value)
    for agent_id, state in new_value.items():
        if agent_id not in log:
            log[agent_id] = []
        log[agent_id].append(state)
    return log


class DeathLogger(Deriver):
    def ports_schema(self):
        return {
            'source': {
                '*': {
                    'boundary': {
                        'death': {}
                    }
                }
            },
            'log': {
                '_default': {},
                '_updater': append_log,
            }}

    def next_update(self, timestep, states):
        source = states['source']
        log = {}
        for agent_id, state in source.items():
            log[agent_id] = state['boundary']['death']
        return {'log': log}
