"""
============
Death Logger
============

This composite can be merged with an environmental simulation in order to save
the final time and state of all agents upon their death.
"""

from vivarium.core.composer import Composer
from vivarium.core.process import Deriver
from vivarium.library.dict_utils import deep_merge
from vivarium.processes.clock import Clock


class DeathLogger(Composer):
    defaults = {'time_step': 1}
    def __init__(self, config=None):
        super().__init__(config)
    def generate_processes(self, config):
        return {
            'clock': Clock({
                'time_step': config['time_step']}),
            'death_log': Logger()}
    def generate_topology(self, config):
        return {
            'clock': {
                'global_time': ('global_time',),
            },
            'death_log': {
                'time': ('global_time',),
                'source': ('agents',),
                'log': ('log',),
            }}


def append_log(current_value, new_value):
    log = deep_merge(dict(current_value), new_value)
    return log


class Logger(Deriver):
    """ Saves the most recent death state."""
    def ports_schema(self):
        return {
            'time': {
                '_default': 0.0,
            },
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
                '_emit': True,
            }}

    def next_update(self, timestep, states):
        source = states['source']
        time = states['time']
        log = {
            agent_id: (time, state['boundary']['death'])
            for agent_id, state in source.items()}
        return {'log': log}
