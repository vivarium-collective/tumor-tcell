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
import math

# vivarium-core imports
from vivarium.core.engine import Engine, timestamp
from vivarium.library.units import units, remove_units
from vivarium.core.control import Control

# plots
from vivarium.plots.agents_multigen import plot_agents_multigen
from tumor_tcell.plots.video import make_video
from tumor_tcell.plots.snapshots import plot_snapshots, format_snapshot_data, get_agent_colors

# tumor-tcell imports
from tumor_tcell.composites.tumor_agent import TumorAgent
from tumor_tcell.composites.t_cell_agent import TCellAgent
from tumor_tcell.composites.tumor_microenvironment import TumorMicroEnvironment
from tumor_tcell.composites.death_logger import DeathLogger

# global parameters
PI = math.pi
TIMESTEP = 60
NBINS = [20, 20]
DEPTH = 15  # um
BOUNDS = [200 * units.um, 200 * units.um]

TUMOR_ID = 'tumor'
TCELL_ID = 'tcell'

# parameters for toy experiments
SMALL_BOUNDS = [20*units.um, 20*units.um]
MEDIUM_BOUNDS = [90*units.um, 90*units.um]

# plotting
TAG_COLORS = {
    ('internal', 'cell_state', 'PDL1p'): 'skyblue',
    ('internal', 'cell_state', 'PDL1n'): 'indianred',
    ('internal', 'cell_state', 'PD1p'): 'limegreen',
    ('internal', 'cell_state', 'PD1n'): 'darkorange', }


def get_tcells(number=1, relative_pd1n=0.2, total_pd1n=None):
    if total_pd1n:
        assert isinstance(total_pd1n, int)
        return {
        '{}_{}'.format(TCELL_ID, n): {
            'type': 'tcell',
            'cell_state': 'PD1n' if n < total_pd1n else 'PD1p',
            'TCR_timer': random.uniform(0, 5400),
            'velocity_timer': 0,
            'velocity': 10.0 * units.um/units.min,
            'diameter': 7.5 * units.um,
        } for n in range(number)}
    else:
        assert relative_pd1n <= 1.0
        return {
        '{}_{}'.format(TCELL_ID, n): {
            'type': 'tcell',
            'cell_state': 'PD1n' if random.uniform(0, 1) < relative_pd1n else 'PD1p',
            'TCR_timer': random.uniform(0, 5400),
            'velocity_timer': 0,
            'velocity': 10.0 * units.um/units.min,
            'diameter': 7.5 * units.um,
        } for n in range(number)}


def get_tumors(number=1, relative_pdl1n=0.5):
    return {
        '{}_{}'.format(TUMOR_ID, n): {
            'type': 'tumor',
            'cell_state': 'PDL1n' if random.uniform(0, 1) < relative_pdl1n else 'PDL1p',
            'diameter': 15 * units.um,
        } for n in range(number)}


def random_location(
        bounds,
        center=None,
        distance_from_center=None,
        excluded_distance_from_center=None
):
    if distance_from_center and excluded_distance_from_center:
        assert distance_from_center > excluded_distance_from_center, \
            'distance_from_center must be greater than excluded_distance_from_center'

    # get the center
    if center:
        center_x = center[0]
        center_y = center[1]
    else:
        center_x = bounds[0]/2
        center_y = bounds[1]/2

    if distance_from_center:
        if excluded_distance_from_center:
            ring_size = distance_from_center - excluded_distance_from_center
            distance = excluded_distance_from_center + ring_size * math.sqrt(random.random())
        else:
            distance = distance_from_center * math.sqrt(random.random())

        angle = random.uniform(0, 2 * PI)
        dy = math.sin(angle)*distance
        dx = math.cos(angle)*distance
        pos_x = center_x+dx
        pos_y = center_y+dy

    elif excluded_distance_from_center:
        in_center = True
        while in_center:
            pos_x = random.uniform(0, bounds[0])
            pos_y = random.uniform(0, bounds[1])
            distance = (pos_x**2 + pos_y**2)**0.5
            if distance > excluded_distance_from_center:
                in_center = False
    else:
        pos_x = random.uniform(0, bounds[0])
        pos_y = random.uniform(0, bounds[1])

    return [pos_x, pos_y]



def lymph_node_location(bounds, location=[[0.95,1],[0.95,1]]):
    return [
        random.uniform(bounds[0]*location[0][0], bounds[0]*location[0][1]),
        random.uniform(bounds[0]*location[1][0], bounds[0]*location[1][1])]

def convert_to_hours(data):
    """Convert seconds to hours"""
    times = list(data.keys())
    for time in times:
        hour = time/3600
        data[hour] = data.pop(time)
    return data


# make defaults
N_TUMORS = 120
N_TCELLS = 9
DEFAULT_TUMORS = get_tumors(number=N_TUMORS)
DEFAULT_TCELLS = get_tcells(number=N_TCELLS)


# The main simulation function
def tumor_tcell_abm(
    bounds=BOUNDS,
    n_bins=NBINS,
    depth=DEPTH,
    field_molecules=['IFNg'],
    n_tumors=120,
    n_tcells=9,
    tumors=None,
    tcells=None,
    tumors_state_PDL1n=0.5,
    tcells_state_PD1n=0.8,
    tcells_total_PD1n=None,
    total_time=70000,
    sim_step=10*TIMESTEP,  # simulation increments at which halt_threshold is checked
    halt_threshold=300,  # stop simulation at this number
    time_step=TIMESTEP,
    emit_step=None,
    emitter='timeseries',
    parallel=False,
    tumors_distance=None,
    tcells_distance=None,
    tumors_excluded_distance=None,
    tcells_excluded_distance=None,
    tumors_center=None,
    tcell_center=None,
    lymph_nodes=False,
):
    """ Tumor-Tcell simulation function

    :param bounds:
    :param n_bins:
    :param depth:
    :param field_molecules:
    :param n_tumors:
    :param n_tcells:
    :param tumors:
    :param tcells:
    :param tumors_state_PDL1n:
    :param tcells_state_PD1n:
    :param total_time:
    :param sim_step:
    :param halt_threshold:
    :param time_step:
    :param emit_step:
    :param emitter:
    :param parallel:
    :param tumors_distance:
    :param tcells_distance:
    :return:
        Simulation output data (dict)
    """
    initial_env_config = {
        'diffusion_field': {'uniform': 0.0}}
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

    #Make initial cells
    if not tcells:
        tcells = get_tcells(
            number=n_tcells,
            relative_pd1n=tcells_state_PD1n,
            total_pd1n=tcells_total_PD1n)

    if not tumors:
        tumors = get_tumors(
            number=n_tumors,
            relative_pdl1n=tumors_state_PDL1n)

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
                    bounds,
                    center=tcell_center,
                    distance_from_center=tcells_distance,
                    excluded_distance_from_center=tcells_excluded_distance,
                ) if not lymph_nodes else lymph_node_location(bounds)),
                'diameter': state.get('diameter', 7.5 * units.um),
                'velocity': state.get('velocity', 10.0 * units.um/units.min)},
            'globals': {
                'LN_no_migration': lymph_nodes,
            },
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
                    bounds,
                    center=tumors_center,
                    distance_from_center=tumors_distance,
                    excluded_distance_from_center=tumors_excluded_distance)),
                'diameter': state.get('diameter', 15 * units.um),
                'velocity': state.get('velocity', 0.0 * units.um/units.min)},
            'neighbors': {
                'present': {
                    'PDL1': state.get('PDL1', None),
                    'MHCI': state.get('MHCI', 1000)}
            }} for agent_id, state in tumors.items()}

    initial_state = {
        **initial_env,
        'agents': {
            **initial_t_cells, **initial_tumors}}

    # make the experiment
    experiment_id = (f"tumor_tcell_{timestamp()}")
    experiment_config = {
        'description': f"n_tcells: {n_tcells} \n"
                       f"n_tumors: {n_tumors} \n"
                       f"tumors_state_PDL1n: {tumors_state_PDL1n} \n"
                       f"tcells_state_PD1n:{tcells_state_PD1n} \n"
                       f"tcells_total_PD1n:{tcells_total_PD1n} \n"
                       f"total_time:{total_time} \n"
                       f"time_step:{time_step} \n"
                       f"sim_step:{sim_step} \n"
                       f"bounds:{bounds} \n"
                       f"n_bins:{n_bins} \n"
                       f"halt_threshold:{halt_threshold} \n"
                       f"tumors_distance:{tumors_distance} \n"
                       f"tcells_distance: {tcells_distance} \n",
        'processes': composite_model.processes,
        'topology': composite_model.topology,
        'initial_state': initial_state,
        'display_info': False,
        'experiment_id': experiment_id,
        'emit_step': emit_step,
        'emitter': {'type': emitter}}
    print(f'Initializing experiment {experiment_id}')
    experiment = Engine(experiment_config)

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
    data = convert_to_hours(data)
    return data


FULL_BOUNDS = [1200*units.um, 1200*units.um] #Usually 1200 by 1200
def full_experiment(
        n_tcells=12,
        n_tumors=1200,
        tcells_state_PD1n=None,
        tumors_state_PDL1n=None,
        tcells_total_PD1n=None,
        lymph_nodes=False,
        total_time=259200,
):

    return tumor_tcell_abm(
        n_tcells=n_tcells,
        n_tumors=n_tumors,
        tcells_state_PD1n=tcells_state_PD1n,
        tumors_state_PDL1n=tumors_state_PDL1n,
        tcells_total_PD1n=tcells_total_PD1n,
        total_time=total_time,
        time_step=TIMESTEP,
        sim_step=100*TIMESTEP,
        emit_step=10*TIMESTEP,
        bounds=FULL_BOUNDS,
        n_bins=[120, 120], #10 um bin size, usually 120 by 120
        halt_threshold=5000, #sqrt(halt_threshold)*15 <bounds, normally 4800
        emitter='database',
        tumors_distance=260*units.um, #sqrt(n_tumors)*15(diameter)/2
        tcells_distance=None, #200*units.um, #in or out (None) of the tumor
        #parallel=True,
        lymph_nodes=lymph_nodes,
    )

#Change experimental PD1 and PDL1 levels for full experiment
def full_experiment_2():
    return full_experiment(
        n_tcells=12,
        tcells_state_PD1n=0.8, #0.2 and 0.8
        tumors_state_PDL1n=0.5, #0.5 originally
        tcells_total_PD1n=9,
    )

def lymph_node_experiment():
    return full_experiment(
        n_tcells=12,
        n_tumors=120,
        tcells_state_PD1n=0.8, #0.2 and 0.8
        tumors_state_PDL1n=0.5, #0.5 originally
        tcells_total_PD1n=9,
        lymph_nodes=True,
        total_time=600000,
    )



def plots_suite(
        data,
        out_dir=None,
        bounds=BOUNDS,
        n_snapshots=8,
        final_time=None,
):

    # separate out tcell and tumor data for multigen plots
    tcell_data = {}
    tumor_data = {}
    for time, time_data in data.items():
        if 'agents' not in time_data:
            continue
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
        'time_display': '(hr)',
        'skip_paths': [
            ('boundary', 'diameter'),
            ('boundary', 'location')]}
    fig1 = plot_agents_multigen(tcell_data, plot_settings, out_dir, TCELL_ID)
    fig2 = plot_agents_multigen(tumor_data, plot_settings, out_dir, TUMOR_ID)

    # snapshots plot
    # extract data
    agents, fields = format_snapshot_data(data)

    # set tag colors.

    fig3 = plot_snapshots(
        bounds=remove_units(bounds),
        agents=remove_units(agents),
        fields=fields,
        tag_colors=TAG_COLORS,
        n_snapshots=n_snapshots,
        final_time=final_time,
        out_dir=out_dir,
        filename='snapshots.pdf',
        default_font_size=48,
        field_label_size=48,
        time_display='hr')

    return fig1, fig2, fig3


def make_snapshot_video(
        data,
        bounds,
        #step=1,   # make frame every n saved steps
        n_steps=10,
        out_dir=None
):
    n_times = len(data.keys())
    step = math.ceil(n_times/n_steps)

    make_video(
        data=remove_units(data),
        bounds=remove_units(bounds),
        agent_shape='circle',
        tag_colors=TAG_COLORS,
        step=step,
        out_dir=out_dir,
        time_display='hr',
        filename='tumor_tcell_video'
    )



# libraries of experiments and plots for easy access by Control
experiments_library = {
    '1': tumor_tcell_abm,
    '4': full_experiment,
    '5': full_experiment_2,
    '6': lymph_node_experiment,
}
plots_library = {
    '1': plots_suite,
    'video': make_snapshot_video,
}
workflow_library = {
    '1': {
        'name': 'tumor_tcell_experiment',
        'experiment': {
            'experiment_id': '1',
            'total_time': 40000,
        },
        'plots': [
            {
                'plot_id': '1',
                'bounds': BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': BOUNDS,
                'n_steps': 100,
            },
        ],
    },
    '2': {
        'name': 'small_experiment',
        'experiment': {
            'experiment_id': '1',
            'bounds': SMALL_BOUNDS,
            'n_tumors': 1,
            'n_tcells': 1,
            'total_time': 10000,
            'n_bins': [1, 1]
        },
        'plots': [
            {
                'plot_id': '1',
                'bounds': SMALL_BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': SMALL_BOUNDS,
            },
        ],
    },
    '3': {
        'name': 'medium_experiment',
        'experiment': {
            'experiment_id': '1',
            'bounds': MEDIUM_BOUNDS,
            'n_tumors': 3,
            'n_tcells': 3,
            'total_time': 50000,
            'n_bins': [3, 3],
            # 'emitter': 'database',
            'tumors_distance': 25 * units.um,
            'tcells_distance': 10 * units.um,
            'parallel': True,
        },
        'plots': [
            {
                'plot_id': '1',
                'bounds': MEDIUM_BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': MEDIUM_BOUNDS,
                'n_steps': 100,
            },
        ],
    },
    '4': {
        'name': 'full_experiment',
        'experiment': '4',
        'plots': [
            {
                'plot_id': '1',
                'bounds': FULL_BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': FULL_BOUNDS,
                'n_steps': 100
            },
        ],
    },
    '5': {
        'name': 'full_experiment_2',
        'experiment': '5',
        'plots': [
            {
                'plot_id': '1',
                'bounds': FULL_BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': FULL_BOUNDS,
                'n_steps': 100
            },
        ],
    },
    '6': {
        'name': 'lymph_node_experiment',
        'experiment': '6',
        'plots': [
            {
                'plot_id': '1',
                'bounds': FULL_BOUNDS
            },
            {
                'plot_id': 'video',
                'bounds': FULL_BOUNDS,
                'n_steps': 100,
            },
        ],
    },
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
    )
