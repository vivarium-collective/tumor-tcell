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
import time as clock
from tqdm import tqdm

# vivarium-core imports
from vivarium.core.experiment import Experiment, timestamp
from vivarium.library.units import units, remove_units
from vivarium.core.control import Control

# plots
from vivarium.plots.agents_multigen import plot_agents_multigen
from tumor_tcell.plots.snapshots import plot_snapshots, format_snapshot_data, get_agent_colors

# tumor-tcell imports
from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment
from tumor_tcell.composites.death_logger import DeathLogger

# global parameters
TIMESTEP = 60
NBINS = [20, 20]
DEPTH = 15  # um
BOUNDS = [200 * units.um, 200 * units.um]

TUMOR_ID = 'tumor'
TCELL_ID = 'tcell'


def get_tcells(number=1, state_per=0.8):
    return {
    '{}_{}'.format(TCELL_ID, n): {
        'type': 'tcell',
        'cell_state': 'PD1n' if random.uniform(0, 1) < state_per else 'PD1p',
        'TCR_timer': random.uniform(0, 5400),
        'velocity_timer': random.uniform(0, 240),
        'velocity': 10.0 * units.um/units.min,
        'diameter': 7.5 * units.um,
    } for n in range(number)}


def get_tumors(number=1, state_per=0.5):
    return {
        '{}_{}'.format(TUMOR_ID, n): {
            'type': 'tumor',
            'cell_state': 'PDL1n' if random.uniform(0, 1) < state_per else 'PDL1p',
            'diameter': 15 * units.um,
        } for n in range(number)}


def random_location(bounds, distance_from_center=None):
    if distance_from_center:
        center_x = bounds[0]/2
        center_y = bounds[1]/2
        return [
            random.uniform(center_x-distance_from_center, center_x+distance_from_center),
            random.uniform(center_y-distance_from_center, center_y+distance_from_center)]
    return [
        random.uniform(0, bounds[0]),
        random.uniform(0, bounds[1])]


# make defaults
N_TUMORS = 120
N_TCELLS = 9
DEFAULT_TUMORS = get_tumors(number=N_TUMORS)
DEFAULT_TCELLS = get_tcells(number=N_TCELLS)


# simulation # 1
def tumor_tcell_abm(
    bounds=BOUNDS,
    n_bins=NBINS,
    depth=DEPTH,
    field_molecules=['IFNg'],
    tumors=DEFAULT_TUMORS,
    tcells=DEFAULT_TCELLS,
    total_time=70000,
    sim_step=10*TIMESTEP,  # simulation increments at which halt_threshold is checked
    halt_threshold=300,  # stop simulation at this number
    time_step=TIMESTEP,
    emit_step=None,
    emitter='timeseries',
    parallel=False,
    tumors_distance=None,
    tcell_distance=None,
):
    initial_env_config = {'uniform': 0.0}
    jitter_force = 0
    t_cell_config = {
        'tcell': {'_parallel': parallel},
        'time_step': time_step}
    tumor_config = {
        'tumor': {'_parallel': parallel},
        'time_step': time_step}
    environment_config = {
        'neighbors_multibody': {
            '_parallel': parallel,
            'time_step': time_step,
            'bounds': bounds,
            'jitter_force': jitter_force},
        'diffusion_field': {
            '_parallel': parallel,
            'time_step': time_step,
            'molecules': field_molecules,
            'bounds': bounds,
            'n_bins': n_bins,
            'depth': depth}}
    logger_config = {'time_step': time_step}

    # make the composers
    t_cell_composer = TCellAgent(t_cell_config)
    tumor_composer = TumorAgent(tumor_config)
    environment_composer = TumorMicroEnvironment(environment_config)
    logger_composer = DeathLogger(logger_config)

    # initialize the composite, and add the environment
    composite_model = logger_composer.generate()
    environment = environment_composer.generate()
    composite_model.merge(composite=environment)

    # add tcells to the composite
    for agent_id in tcells.keys():
        t_cell = t_cell_composer.generate({'agent_id': agent_id})
        composite_model.merge(composite=t_cell, path=('agents', agent_id))

    # add tumors to the composite
    for agent_id in tumors.keys():
        tumor = tumor_composer.generate({'agent_id': agent_id})
        composite_model.merge(composite=tumor, path=('agents', agent_id))

    # make initial environment state
    initial_env = composite_model.initial_state(initial_env_config)

    # initialize cell state
    initial_t_cells = {
        agent_id: {
            'boundary': {
                'location': state.get('location', random_location(
                    bounds, distance_from_center=tcell_distance)),
                'diameter': state.get('diameter', 10 * units.um),
                'velocity': state.get('velocity', 0.0 * units.um/units.min)},
            'internal': {
                'cell_state': state.get('cell_state', None),
                'velocity_timer': state.get('velocity_timer', 0),
                'TCR_timer': state.get('TCR_timer', 0)},
            'neighbors': {
                'present': {
                    'PD1': state.get('PD1', None),
                    'TCR': state.get('TCR', 50000)}
            }} for agent_id, state in tcells.items()}

    initial_tumors = {
        agent_id: {
            'internal': {
                'cell_state': state.get('cell_state', None)},
            'boundary': {
                'location': state.get('location', random_location(
                    bounds, distance_from_center=tumors_distance)),
                'diameter': state.get('diameter', 20 * units.um),
                'velocity': state.get('velocity', 0.0 * units.um/units.min)},
            'neighbors': {
                'present': {
                    'PDL1': state.get('PDL1', None),
                    'MHCI': state.get('MHCI', None)}
            }} for agent_id, state in tumors.items()}

    initial_state = {
        **initial_env,
        'agents': {
            **initial_t_cells, **initial_tumors}}

    # make the experiment
    experiment_id = (f"tumor_tcell_{timestamp()}")
    experiment_config = {
        'processes': composite_model.processes,
        'topology': composite_model.topology,
        'initial_state': initial_state,
        'display_info': False,
        'experiment_id': experiment_id,
        'emit_step': emit_step,
        'emitter': {'type': emitter}}
    print(f'Initializing experiment {experiment_id}')
    experiment = Experiment(experiment_config)

    # run simulation and terminate upon reaching total_time or halt_threshold
    clock_start = clock.time()
    for time in tqdm(range(0, total_time, sim_step)):
        n_agents = len(experiment.state.get_value()['agents'])
        if n_agents < halt_threshold:
            experiment.update(sim_step)
        else:
            print(f'halt threshold of {halt_threshold} reached at time = {time}')
            continue

    # print runtime and finalize
    clock_finish = clock.time() - clock_start
    print('Completed in {:.2f} seconds'.format(clock_finish))
    experiment.end()

    # return time
    data = experiment.emitter.get_data_deserialized()
    return data


FULL_BOUNDS = [1200*units.um, 1200*units.um]
def full_experiment():
    return tumor_tcell_abm(
        # tumors=get_tumors(number=3500),
        # tcells=get_tcells(number=30),
        # total_time=259200,
        tumors=get_tumors(number=1200),
        tcells=get_tcells(number=12),
        total_time=259200,
        time_step=TIMESTEP,
        sim_step=100*TIMESTEP,
        emit_step=10*TIMESTEP,
        bounds=FULL_BOUNDS,
        n_bins=[120, 120], #10 um bin size
        halt_threshold=4800, #sqrt(halt_threshold)*15 <bounds
        emitter='database',
        tumors_distance=260*units.um, #sqrt(n_tumors)*15(diameter)/2
        tcell_distance=200*units.um, #in or out of the tumor
        #parallel=True,
    )


MEDIUM_BOUNDS = [90*units.um, 90*units.um]
def medium_experiment():
    return tumor_tcell_abm(
        tumors=get_tumors(number=3),
        tcells=get_tcells(number=3),
        total_time=50000,
        bounds=MEDIUM_BOUNDS,
        n_bins=[3, 3],
        emitter='database',
        tumors_distance=25 * units.um,
        tcell_distance=10 * units.um,
    )


SMALL_BOUNDS = [20*units.um, 20*units.um]
def small_experiment():
    return tumor_tcell_abm(
        tumors=get_tumors(number=1),
        tcells=get_tcells(number=1),
        total_time=100000,
        bounds=SMALL_BOUNDS,
        n_bins=[1, 1])


def plots_suite_small_bounds(
        data, out_dir=None, bounds=SMALL_BOUNDS):
    return plots_suite(data, out_dir, bounds)


def plots_suite_medium_bounds(
        data, out_dir=None, bounds=MEDIUM_BOUNDS):
    return plots_suite(data, out_dir, bounds)


def plots_suite_full_bounds(
        data, out_dir=None, bounds=FULL_BOUNDS):
    return plots_suite(data, out_dir, bounds)


def plots_suite(
        data, out_dir=None, bounds=BOUNDS):

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

    # get the final death log
    times_vector = list(data.keys())
    death_log = data[times_vector[-1]]['log']

    # make multigen plot for tcells and tumors
    plot_settings = {
        'skip_paths': [
            ('boundary', 'diameter'),
            ('boundary', 'location')]}
    fig1 = plot_agents_multigen(tcell_data, plot_settings, out_dir, TCELL_ID)
    fig2 = plot_agents_multigen(tumor_data, plot_settings, out_dir, TUMOR_ID)

    # snapshots plot
    # extract data
    agents, fields = format_snapshot_data(data)

    # set tag colors.
    tag_colors = {
        ('internal', 'cell_state', 'PDL1p'): 'skyblue',
        ('internal', 'cell_state', 'PDL1n'): 'indianred',
        ('internal', 'cell_state', 'PD1p'): 'limegreen',
        ('internal', 'cell_state', 'PD1n'): 'darkorange',}

    fig3 = plot_snapshots(
        bounds=remove_units(bounds),
        agents=remove_units(agents),
        fields=fields,
        tag_colors=tag_colors,
        n_snapshots=8,
        out_dir=out_dir,
        filename='snapshots')

    return fig1, fig2, fig3


# libraries of experiments and plots for easy access by Control
experiments_library = {
    '1': tumor_tcell_abm,
    '2': small_experiment,
    '3': medium_experiment,
    '4': full_experiment,
}
plots_library = {
    '1': plots_suite,
    '2': plots_suite_small_bounds,
    '3': plots_suite_medium_bounds,
    '4': plots_suite_full_bounds,
}
workflow_library = {
    '1': {
        'name': 'tumor_tcell_experiment',
        'experiment': '1',
        'plots': ['1'],
    },
    '2': {
        'name': 'small_experiment',
        'experiment': '2',
        'plots': ['2'],
    },
    '3': {
        'name': 'medium_experiment',
        'experiment': '3',
        'plots': ['3'],
    },
    '4': {
        'name': 'full_experiment',
        'experiment': '4',
        'plots': ['4'],
    },
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
    )
