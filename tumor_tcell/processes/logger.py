from vivarium.core.process import Deriver

def append_update(current_value, new_value):
    assert isinstance(current_value, list)
    return current_value.append(new_value)

class Logger(Deriver):
    defaults = {
        'source_schema': {}}
    def ports_schema(self):
        return {
            'source': self.parameters['source_schema'],
            'log': {
                '_default': [],
                '_updater': append_update,
            }}
    def next_update(self, timestep, states):
        state = states['source']

        import ipdb; ipdb.set_trace()

        return {'log': state}
