"""
========================
Tumor/T-cell Experiments
========================

Experiments can be triggered from the command line:

```
$ python tumor_tcell/experiments/main.py [experiment_name]
```
"""

import random

# vivarium-core imports
from vivarium.core.composition import (
    agent_environment_experiment,
    make_agent_ids,
    compose_experiment,
    GENERATORS_KEY,
    EXPERIMENT_OUT_DIR,
)
from vivarium.library.units import units, remove_units
from vivarium.core.control import Control

# plots
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_cell.plots.multibody_physics import plot_tags, plot_snapshots

# tumor-tcell imports
from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment


out_dir = EXPERIMENT_OUT_DIR
TUMOR_ID = 'tumor'
TCELL_ID = 'tcell'

# TODO -- pass this in through config
BOUNDS = [200 * units.um, 200 * units.um]


# simulation # 1
def simulation_1():
    # experiment parameters
    total_time = 10000
    bounds = BOUNDS
    time_step = 60

    # cell config
    n_tcells = 5
    n_tumors = 5
    t_cell_ids = [TUMOR_ID + '_' + str(num) for num in range(n_tcells)]
    tumor_ids = [TCELL_ID + '_' + str(num) for num in range(n_tumors)]

    # initial state
    initial_t_cells = {
        agent_id: {
            'boundary': {
                'location': [
                    random.uniform(0, bounds[0]),
                    random.uniform(0, bounds[1])]
            }
        } for agent_id in t_cell_ids
    }
    initial_tumors = {
        agent_id: {
            'boundary': {
                'location': [
                    random.uniform(0, bounds[0]),
                    random.uniform(0, bounds[1])]
            }
        } for agent_id in tumor_ids
    }
    initial_state = {
        'fields': {},
        'agents': {
            **initial_t_cells, **initial_tumors}
    }

    # t-cell configurations
    t_cell_hierarchy = {
        agent_id: {
            GENERATORS_KEY: {
                'type': TCellAgent,
                'config': {
                    'time_step': time_step,
                    'agent_id': agent_id,
                }
            }
        } for agent_id in t_cell_ids
    }

    # tumor configurations
    tumor_hierarchy = {
        agent_id: {
            GENERATORS_KEY: {
                'type': TumorAgent,
                'config': {
                    'time_step': time_step,
                    'agent_id': agent_id,
                }
            }
        } for agent_id in tumor_ids
    }

    # declare the hierarchy full
    hierarchy = {
        # generate the tumor micro-environment at the top level
        GENERATORS_KEY: {
            'type': TumorMicroEnvironment,
            'config': {
                'neighbors_multibody': {
                    'time_step': time_step,
                    'bounds': bounds,
                },
                'diffusion_field': {
                    'time_step': time_step,
                    'bounds': bounds,
                }
            }
        },
        # cells are one level down, under the 'agents' key
        'agents': {**t_cell_hierarchy, **tumor_hierarchy}
    }

    # configure experiment
    experiment = compose_experiment(
        hierarchy=hierarchy,
        initial_state=initial_state
    )

    # run simulation
    experiment.update(total_time)
    data = experiment.emitter.get_data()
    experiment.end()

    return data


def plots_suite(data, out_dir=EXPERIMENT_OUT_DIR):

    # separate out tcell and tumor data for multigen plots
    tcell_data = {}
    tumor_data = {}
    for time, time_data in data.items():
        all_agents_data = time_data['agents']
        tcell_data[time] = {
            'agents': {
                agent_id: agent_data
                for agent_id, agent_data in all_agents_data.items()
                if TCELL_ID in agent_id}}
        tumor_data[time] = {
            'agents': {
                agent_id: agent_data
                for agent_id, agent_data in all_agents_data.items()
                if TUMOR_ID in agent_id}}

    # make multigen plot for tcells and tumors
    plot_settings = {}
    plot_agents_multigen(tcell_data, plot_settings, out_dir, TCELL_ID)
    plot_agents_multigen(tumor_data, plot_settings, out_dir, TUMOR_ID)

    # snapshots plot
    # extract data
    env_config = {'bounds': remove_units(BOUNDS)}
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    # fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_data = {
        'agents': agents,
        'fields': {},
        'config': env_config}
    plot_config = {
        'fields': [],
        'n_snapshots': 8,
        'agent_shape': 'circle',
        'out_dir': out_dir}
    plot_snapshots(plot_data, plot_config)



# all of the experiments go here for easy access by control class
experiments_library = {
    '1': simulation_1,
}
plots_library = {
    '1': plots_suite,
}
workflow_library = {
    '1': {
        'name': 'tumor_tcell_experiment',
        'experiment': '1',
        'plots': ['1'],
    }
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
    )