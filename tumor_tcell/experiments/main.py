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

# global parameters
BOUNDS = [100 * units.um, 100 * units.um]

TUMOR_ID = 'tumor'
TCELL_ID = 'tcell'
N_TUMORS = 14
N_TCELLS = 9
DEFAULT_TUMORS = {
    '{}_{}'.format(TUMOR_ID, n): {
        #'location': [x, y],
        'type': 'tumor',
        'cell_state': 'PDL1n' if random.uniform(0, 1) < 0.5 else 'PDL1p',
        #'PDL1': 50000,
        #'MHCI': 50000,
        'diameter': 20,
    } for n in range(N_TUMORS)
}
DEFAULT_TCELLS = {
    '{}_{}'.format(TCELL_ID, n): {
        #'location': [x, y],
        'type': 'tcell',
        'cell_state': 'PD1n',
        #'PD1': 0,
        #'TCR':50000,
        'diameter': 10,
    } for n in range(N_TCELLS)
}

def random_location(bounds):
    return [
        random.uniform(0, bounds[0]),
        random.uniform(0, bounds[1])]


# simulation # 1
def tumor_tcell_abm(
    bounds=BOUNDS,
    field_molecules=['IFNg'],
    tumors=DEFAULT_TUMORS,
    tcells=DEFAULT_TCELLS,
    total_time=2000,
    time_step=60
):

    ## configure the cells
    # t-cell configuration
    t_cell_hierarchy = {
        agent_id: {
            GENERATORS_KEY: {
                'type': TCellAgent,
                'config': {
                    'time_step': time_step,
                    'agent_id': agent_id,
                }
            }
        } for agent_id in tcells.keys()
    }

    # tumor configuration
    tumor_hierarchy = {
        agent_id: {
            GENERATORS_KEY: {
                'type': TumorAgent,
                'config': {
                    'time_step': time_step,
                    'agent_id': agent_id,
                }
            }
        } for agent_id in tumors.keys()
    }

    # declare the full hierarchy with the environments
    hierarchy = {
        # generate the tumor micro-environment at the top level
        GENERATORS_KEY: {
            'type': TumorMicroEnvironment,
            'config': {
                'neighbors_multibody': {
                    'time_step': time_step,
                    'bounds': bounds,
                    'jitter_force': 5e-4,
                },
                'diffusion_field': {
                    'time_step': time_step,
                    'molecules': field_molecules,
                    'bounds': bounds,
                }
            }
        },
        # cells are one level down, under the 'agents' key
        'agents': {**t_cell_hierarchy, **tumor_hierarchy}
    }

    # make environment instance to get an initial state
    environment = TumorMicroEnvironment(hierarchy[GENERATORS_KEY]['config'])
    initial_env = environment.initial_state()

    # initialize state
    initial_t_cells = {
        agent_id: {
            'boundary': {
                'location': state.get('location', random_location(bounds)),
                'diameter': state.get('diameter', 10) * units.um,
            },
            'internal': {
                'cell_state': state.get('cell_state', None)
            },
            'neighbors': {
                'present': {
                    'PD1': state.get('PD1', None),
                    'TCR': state.get('TCR', 50000)
                }
            },
        } for agent_id, state in tcells.items()
    }
    initial_tumors = {
        agent_id: {
            'internal': {
                'cell_state': state.get('cell_state', None)
            },
            'boundary': {
                'location': state.get('location', random_location(bounds)),
                'diameter': state.get('diameter', 20) * units.um,
            },
            'neighbors': {
                'present': {
                    'PDL1': state.get('PDL1', None),
                    'MHCI': state.get('MHCI', None)
                }
            },
        } for agent_id, state in tumors.items()
    }
    initial_state = {
        **initial_env,
        'agents': {
            **initial_t_cells, **initial_tumors}
    }

    # configure the simulation experiment
    experiment = compose_experiment(
        hierarchy=hierarchy,
        initial_state=initial_state)

    # run simulation
    experiment.update(total_time)
    data = experiment.emitter.get_data()
    experiment.end()

    return data


def plots_suite(data, out_dir=None):

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
    plot_settings = {
        'skip_paths': [
            ('boundary', 'diameter'),
            ('boundary', 'location'),
        ]
    }
    fig1 = plot_agents_multigen(tcell_data, plot_settings, out_dir, TCELL_ID)
    fig2 = plot_agents_multigen(tumor_data, plot_settings, out_dir, TUMOR_ID)

    # snapshots plot
    # extract data
    env_config = {'bounds': remove_units(BOUNDS)}
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_data = {
        'agents': agents,
        'fields': fields,
        'config': env_config}
    plot_config = {
        'fields': [],
        'n_snapshots': 8,
        'agent_shape': 'circle',
        'out_dir': out_dir}
    fig3 = plot_snapshots(plot_data, plot_config)


    return fig1, fig2, fig3

# all of the experiments go here for easy access by control class
experiments_library = {
    '1': tumor_tcell_abm,
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
