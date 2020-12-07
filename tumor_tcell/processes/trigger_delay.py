from vivarium.core.process import Deriver


class DelayTrigger(Deriver):
    def ports_schema(self):
        return {
            'source': {'_default': False},
            'target': {'_default': False}}
    def next_update(self, timestep, states):
        state = states['source']
        if state:
            return {
                'target': {
                    '_updater': 'set',
                    '_value': state}}
        else:
            return {}
