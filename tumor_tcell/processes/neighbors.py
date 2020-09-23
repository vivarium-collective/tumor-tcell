"""
================================================
Multibody physics process with neighbor tracking
================================================
"""

from __future__ import absolute_import, division, print_function

import os
import random
import math

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# vivarium imports
from tumor_tcell.library.pymunk_multibody import PymunkMultibody
from vivarium.library.units import units, remove_units
from vivarium.core.process import Process
from vivarium.core.composition import process_in_experiment

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'neighbors'
DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi



def sphere_volume_from_diameter(diameter):
    radius = diameter / 2
    volume = 4 / 3 * (PI * radius**3)
    return volume

def make_random_position(bounds):
    return [
        np.random.uniform(0, bounds[0]),
        np.random.uniform(0, bounds[1])]


class Neighbors(Process):
    """ Neighbors process for tracking cell bodies.

    Simulates collisions between cell bodies with a physics engine.

    :term:`Ports`:
    * ``cells``: The store containing all cell sub-compartments. Each cell in
      this store has values for location, diameter, mass.

    Arguments:
        parameters(dict): Accepts the following configuration keys:

        * **jitter_force**: force applied to random positions along cell
          bodies to mimic thermal fluctuations. Produces Brownian motion.
        * **bounds** (:py:class:`list`): size of the environment in
          micrometers, with ``[x, y]``.
        * ***animate*** (:py:class:`bool`): interactive matplotlib option to
          animate multibody. To run with animation turned on set True, and use
          the TKAgg matplotlib backend:

          .. code-block:: console

              $ MPLBACKEND=TKAgg python tumor_tcell/processes/neighbors.py
    """

    name = NAME
    defaults = {
        'time_step': 2,
        'cells': {},
        'jitter_force': 0.0,  # pN
        'bounds': DEFAULT_BOUNDS,
        'animate': False,
    }

    def __init__(self, parameters=None):
        super(Neighbors, self).__init__(parameters)

        # multibody parameters
        jitter_force = self.parameters['jitter_force']
        self.bounds = self.parameters['bounds']

        # make the multibody object
        time_step = self.parameters['time_step']
        multibody_config = {
            'cell_shape': 'circle',
            'jitter_force': jitter_force,
            'bounds': self.bounds,
            'physics_dt': time_step / 10,
        }
        self.physics = PymunkMultibody(multibody_config)

        # interactive plot for visualization
        self.animate = self.parameters['animate']
        if self.animate:
            plt.ion()
            self.ax = plt.gca()
            self.ax.set_aspect('equal')

    def ports_schema(self):
        glob_schema = {
            '*': {
                'boundary': {
                    'location': {
                        '_emit': True,
                        '_default': [
                            0.5 * bound for bound in self.bounds],
                        '_updater': 'set',
                        '_divider': 'set'},
                    'diameter': {
                        '_emit': True,
                        '_default': 1.0, #* units.um,
                        '_divider': 'split',
                        '_updater': 'set'},
                    'mass': {
                        '_emit': True,
                        '_default': 1 * units.fg,
                        '_updater': 'set'},
                },
                'neighbors': {
                    '*': {}
                    # 'PD1': {},
                    # 'cytotoxic_packets': {},
                    # 'PDL1': {},
                    # 'MHCI': {}
                    # TODO -- tumors required cytotoxic_packets and PD1
                    # TODO -- t cells require PDL1 and MHCI
                }
            }
        }
        schema = {'cells': glob_schema}
        return schema

    def next_update(self, timestep, states):
        cells = states['cells']

        # animate before update
        if self.animate:
            self.animate_frame(cells)

        # update multibody with new calls
        self.physics.update_bodies(remove_units(cells))

        # run simulation
        self.physics.run(timestep)

        # get new cell positions
        cell_positions = self.physics.get_body_positions()

        # get neighbors
        cell_neighbors = self.get_neighbors(cell_positions)

        # import ipdb; ipdb.set_trace()

        update = {
            'cells': {
                cell_id: {
                    'boundary': {
                        'location': list(cell_positions[cell_id])
                    }
                } for cell_id in cells.keys()
            }
        }
        return update

    def get_neighbors(self, cell_positions):
        cell_neighbors = {}
        for cell_id, position in cell_positions.items():
            other_cell_locations = {
                location: other_id
                for other_id, location in cell_positions.items()
                if other_id is not cell_id}
            dist = lambda x, y: (x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2
            closest = min(list(other_cell_locations.keys()), key=lambda co: dist(co, position))
            neighbor_id = other_cell_locations[closest]
            cell_neighbors[cell_id] = neighbor_id
        return cell_neighbors

    ## matplotlib interactive plot
    def animate_frame(self, cells):
        plt.cla()
        for cell_id, data in cells.items():
            # location, orientation, length
            data = data['boundary']
            x_center = data['location'][0]
            y_center = data['location'][1]
            diameter = data['diameter']

            # get bottom left position
            radius = (diameter / 2)
            x = x_center - radius
            y = y_center - radius

            # Create a circle
            circle = patches.Circle((x, y), radius, linewidth=1, edgecolor='b')
            self.ax.add_patch(circle)

        plt.xlim([0, self.bounds[0]])
        plt.ylim([0, self.bounds[1]])
        plt.draw()
        plt.pause(0.01)


# configs
def single_cell_config(config):
    # cell dimensions
    diameter = 1
    volume = sphere_volume_from_diameter(diameter)
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    location = config.get('location')
    if location:
        location = [loc * bounds[n] for n, loc in enumerate(location)]
    else:
        location = make_random_position(bounds)
    return {
        'boundary': {
        'location': location,
        'volume': volume,
        'diameter': diameter,
        'mass': 1339 * units.fg,
        'thrust': 0,
        'torque': 0}}


def cell_body_config(config):
    cell_ids = config['cell_ids']
    cell_config = {
        cell_id: single_cell_config(config)
        for cell_id in cell_ids}
    return {'cells': cell_config}


# tests and simulations
class InvokeUpdate(object):
    def __init__(self, update):
        self.update = update
    def get(self, timeout=0):
        return self.update


def simulate_growth_division(config, settings):
    initial_cells_state = config['cells']

    # make the process
    multibody = Neighbors(config)
    experiment = process_in_experiment(multibody)
    experiment.state.update_subschema(
        ('cells',), {
            'boundary': {
                'mass': {
                    '_divider': 'split'},
                }})
    experiment.state.apply_subschemas()

    # get initial cell state
    experiment.state.set_value({'cells': initial_cells_state})
    cells_store = experiment.state.get_path(['cells'])

    ## run simulation
    # get simulation settings
    growth_rate = settings.get('growth_rate', 0.0006)
    growth_rate_noise = settings.get('growth_rate_noise', 0.0)
    division_volume = settings.get('division_volume', 0.4)
    total_time = settings.get('total_time', 10)
    timestep = 1

    time = 0
    while time < total_time:
        experiment.update(timestep)
        time += timestep
        cells_state = cells_store.get_value()

        invoked_update = []
        for cell_id, state in cells_state.items():
            state = state['boundary']
            location = state['location']
            diameter = state['diameter']
            mass = state['mass'].magnitude

            # update
            growth_rate2 = (growth_rate + np.random.normal(0.0, growth_rate_noise)) * timestep
            new_mass = mass + mass * growth_rate2
            new_diameter = diameter + diameter * growth_rate2
            new_volume = sphere_volume_from_diameter(new_diameter)

            if new_volume > division_volume:
                daughter_ids = [str(cell_id) + '0', str(cell_id) + '1']
                daughter_updates = []
                for daughter_id in daughter_ids:
                    daughter_updates.append({
                        'daughter': daughter_id,
                        'path': (daughter_id,),
                        'processes': {},
                        'topology': {},
                        'initial_state': {}})
                update = {
                    '_divide': {
                        'mother': cell_id,
                        'daughters': daughter_updates}}
            else:
                update = {
                    cell_id: {
                        'boundary': {
                            'volume': new_volume,
                            'diameter': new_diameter,
                            'mass': new_mass * units.fg}}}

            invoked_update.append((InvokeUpdate({'cells': update}), None, None))

        # update experiment
        experiment.send_updates(invoked_update)

    return experiment.emitter.get_data()

def multibody_neighbors_workflow(config={}, out_dir='out', filename='neighbors'):
    n_cells = 2
    cell_ids = [str(cell_id) for cell_id in range(n_cells)]

    bounds = [20, 20]
    settings = {
        'growth_rate': 0.02,
        'growth_rate_noise': 0.02,
        'division_volume': 2.6,
        'total_time': 120}
    gd_config = {
        'animate': True,
        'jitter_force': 1e0,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'cell_ids': cell_ids}
    gd_config.update(cell_body_config(body_config))
    gd_data = simulate_growth_division(gd_config, settings)



if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    multibody_neighbors_workflow({'animate': True}, out_dir)