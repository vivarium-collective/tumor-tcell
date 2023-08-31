"""
================================================
Multibody physics process with neighbor tracking
================================================

Simulates collisions between cell bodies with a physics engine,
and manages the exchange of a few key molecules.

Each agent is expected to have a 'neighbors' dictionary and a 'boundary' dictionary/
Exchange works by reading the 'neighbors' dictionary for the following keys:
    Typically used for soluble signaling molecules.
    * **transfer**: Giving from neighbor cell through 'receive'.
    * **receive**: Getting from neighbor cell through 'transfer'.
    Typically used for membrane-bound signaling molecules.
    * **accept**: The biological response to the binding action with molecules in 'present'.
    * **present**: The ligand, which triggers a biological response for cells with matching 'accept' molecules.

*Notes*:
    * cell_type must be either 'tumor' or 't_cell'
"""

import os
import random
import math
import copy

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from tumor_tcell.library.location import DEFAULT_BOUNDS, DEFAULT_LENGTH_UNIT
# vivarium imports
from tumor_tcell.library.pymunk_minimal import PymunkMinimal as Pymunk
from vivarium.library.units import units, remove_units
from vivarium.core.process import Process
from vivarium.core.composition import process_in_experiment

# directories
from tumor_tcell import PROCESS_OUT_DIR


NAME = 'neighbors'
DEFAULT_MASS_UNIT = units.ng
DEFAULT_VELOCITY_UNIT = units.um / units.s

# constants
PI = math.pi


# helper functions
def daughter_locations(mother_location, state):
    """Given a mother's location, place daughter cells close"""
    mother_diameter = state['diameter']
    mother_x = mother_location[0]
    mother_y = mother_location[1]
    locations = []
    for daughter in range(2):
        location = [mother_x + random.gauss(0, 0.1)*mother_diameter, mother_y + random.gauss(0, 0.1)*mother_diameter]
        locations.append(location)
    return locations

def sphere_volume_from_diameter(diameter):
    """Given a diameter for a sphere, return the volume"""
    radius = diameter / 2
    volume = 4 / 3 * (PI * radius**3)
    return volume

def make_random_position(bounds):
    """Return a random position with the [x, y] bounds"""
    return [
        np.random.uniform(0, bound.magnitude) * bound.units
        for bound in bounds]

def add_to_dict(dict, added):
    for k, v in added.items():
        if k in dict:
            dict[k] += v
        else:
            dict[k] = v
    return dict

def remove_from_dict(dict, removed):
    for k, v in removed.items():
        if k in dict:
            dict[k] -= v
        else:
            dict[k] = -v
    return dict


class Neighbors(Process):
    """Neighbors process for tracking cell bodies.

    Simulates collisions between cell bodies with a physics engine.

    :term:`Ports`:
    * ``cells``: The store containing all cell sub-compartments. Each cell in
      this store has values for location, diameter, mass.

    Arguments:
        parameters(dict): Accepts the following configuration keys:

        * **jitter** (float): random force applied to cells, in pN.
        * **bounds** (list): size of the environment in the units specified by length_unit, with ``[x, y]``.
        * **length_unit** (Quantity): standard length unit for the physics engine -- this includes the bounds.
        * **mass_unit** (Quantity): standard mass unit for the physics engine.
        * **velocity_unit** (Quantity): standard velocity unit for the physics engine.
        * **neighbor_distance** (Quantity): distance from an agent's outer membrane at which
            other agents are considered neighbors.
        * ***animate*** (:py:class:`bool`): interactive matplotlib option to
          animate multibody. To run with animation turned on set True, and use
          the TKAgg matplotlib backend:

          .. code-block:: console

              $ MPLBACKEND=TKAgg python tumor_tcell/processes/neighbors.py
    """

    name = NAME
    defaults = {
        'jitter_force': 1e-6,
        'bounds': remove_units(DEFAULT_BOUNDS),
        'length_unit': DEFAULT_LENGTH_UNIT,
        'mass_unit': DEFAULT_MASS_UNIT,
        'velocity_unit': DEFAULT_VELOCITY_UNIT,
        'neighbor_distance': 1 * units.um,
        'animate': False,
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        self.length_unit = self.parameters['length_unit']
        self.mass_unit = self.parameters['mass_unit']
        self.velocity_unit = self.parameters['velocity_unit']
        self.neighbor_distance = self.parameters['neighbor_distance'].to(self.length_unit).magnitude
        self.cell_loc_units = {}

        # make the multibody object
        timestep = self.parameters['timestep']
        multibody_config = {
            'cell_shape': 'circle',
            'jitter_force': self.parameters['jitter_force'],
            'bounds': [
                b.to(self.length_unit).magnitude
                for b in parameters['bounds']],
            'physics_dt': min(timestep/10, 0.1)}
        self.physics = Pymunk(multibody_config)

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
                    # cell_type must be either 'tumor' or 't_cell'
                    'cell_type': {},
                    'location': {
                        '_emit': True,
                        '_default': [
                            0.5 * bound for bound in self.parameters['bounds']],
                        '_updater': 'set',
                        '_divider': {
                            'divider': daughter_locations,
                            'topology': {
                                'diameter': ('..', 'diameter',)
                            },
                        }
                    },
                    'diameter': {
                        '_emit': True,
                        '_default': 1.0 * self.length_unit},
                    'mass': {
                        '_default': 1.0 * self.mass_unit},
                    'velocity': {
                        '_default': 0.0 * self.velocity_unit,
                    }
                },
                'neighbors': {
                    'present': {
                        '*': {
                            '_default': 0.0,
                            '_updater': 'set',
                        }
                    },
                    'accept': {
                        '*': {
                            '_default': 0.0,
                            '_updater': 'set',
                        }
                    },
                    'transfer': {
                        '*': {
                            '_default': 0.0,
                        }
                    },
                    'receive': {
                        '*': {
                            '_default': 0.0,
                        }
                    },
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

        # update physics with new cells
        cells = self.bodies_remove_units(cells)
        self.physics.update_bodies(cells)

        # run simulation
        self.physics.run(timestep)

        # get new cell positions and neighbors
        cell_positions = self.physics.get_body_positions()
        cell_neighbors = self.get_all_neighbors(cells, cell_positions)

        # add units to cell_positions
        cell_positions = self.location_add_units(cell_positions)

        # exchange with neighbors
        exchange = {
            cell_id: {
                'accept': {},
                'present': {},
                'transfer': {},
                'receive': {},
            } for cell_id in cells.keys()}

        for cell_id, neighbors in cell_neighbors.items():
            for neighbor_id in neighbors:
                # the neighbor's transfer moves to the cell's receive and then is removed
                transfer = cells[cell_id]['neighbors']['transfer']
                exchange[neighbor_id]['receive'] = add_to_dict(exchange[neighbor_id]['receive'], transfer)
                exchange[cell_id]['transfer'] = remove_from_dict(exchange[cell_id]['transfer'], transfer)

                # present and accept are not removed but updated for each other
                present = cells[neighbor_id]['neighbors']['present']
                exchange[cell_id]['accept'] = add_to_dict(exchange[cell_id]['accept'], present)

        update = {
            'cells': {
                cell_id: {
                    'boundary': {
                        'location': list(cell_positions[cell_id])},
                    'neighbors': exchange[cell_id]
                } for cell_id in cells.keys()
            }
        }

        return update

    def bodies_remove_units(self, bodies):
        """Convert units to the standard, and then remove them

        This is required for interfacing the physics engine, which does not track units
        """
        for bodies_id, specs in bodies.items():
            # convert location
            bodies[bodies_id]['boundary']['location'] = [loc.to(self.length_unit).magnitude for loc in specs['boundary']['location']]
            # convert diameter
            bodies[bodies_id]['boundary']['diameter'] = specs['boundary']['diameter'].to(self.length_unit).magnitude
            # convert mass
            bodies[bodies_id]['boundary']['mass'] = specs['boundary']['mass'].to(self.mass_unit).magnitude
            # convert velocity
            bodies[bodies_id]['boundary']['velocity'] = specs['boundary']['velocity'].to(self.velocity_unit).magnitude
        return bodies

    def location_add_units(self, bodies):
        for body_id, location in bodies.items():
            bodies[body_id] = [(loc * self.length_unit) for loc in location]
        return bodies

    def get_neighbors(self, cell_loc, cell_radius, neighbor_loc, neighbor_radius):
        neighbors = {}
        for neighbor_id, loc in neighbor_loc.items():
            distance = ((cell_loc[0] - loc[0]) ** 2 + (cell_loc[1] - loc[1]) ** 2) ** 0.5
            neighbor_rad = neighbor_radius[neighbor_id]
            inner_distance = distance - cell_radius - neighbor_rad
            if inner_distance <= self.neighbor_distance:
                neighbors[neighbor_id] = inner_distance
        return neighbors

    def get_all_neighbors(self, cells, current_positions):
        """
        only count neighbor if they are within 'neighbor_distance' from outer boundary of cell
        """

        tcell_positions = {
            cell_id: current_positions[cell_id]
            for cell_id, specs in cells.items()
            if specs['boundary']['cell_type'] == 't-cell'}
        tumor_positions = {
            cell_id: current_positions[cell_id]
            for cell_id, specs in cells.items()
            if specs['boundary']['cell_type'] == 'tumor'}
        cell_radii = {
            cell_id: (specs['boundary']['diameter'] / 2)
            for cell_id, specs in cells.items()}

        cell_neighbors = {}

        # t-cells polarize to one tumor cell: find the closest
        for cell_id, location in tcell_positions.items():
            radius = cell_radii[cell_id]
            neighbors = self.get_neighbors(location, radius, tumor_positions, cell_radii)
            if neighbors:
                cell_neighbors[cell_id] = [min(neighbors, key=neighbors.get)]
            else:
                cell_neighbors[cell_id] = []

        # tumors can have multiple t-cell neighbors
        for cell_id, location in tumor_positions.items():
            radius = cell_radii[cell_id]
            neighbors = self.get_neighbors(location, radius, tcell_positions, cell_radii)
            if neighbors:
                cell_neighbors[cell_id] = list(neighbors.keys())
            else:
                cell_neighbors[cell_id] = []

        return cell_neighbors

    def remove_length_units(self, value):
        return value.to(self.length_unit).magnitude

    def animate_frame(self, cells):
        """matplotlib interactive plot"""
        plt.cla()
        bounds = copy.deepcopy(self.parameters['bounds'])
        for cell_id, data in cells.items():
            # location, orientation, length
            data = data['boundary']
            x_center = self.remove_length_units(data['location'][0])
            y_center = self.remove_length_units(data['location'][1])
            diameter = self.remove_length_units(data['diameter'])

            # get bottom left position
            radius = (diameter / 2)
            x = x_center - radius
            y = y_center - radius

            # Create a circle
            circle = patches.Circle((x, y), radius, linewidth=1, edgecolor='b')
            self.ax.add_patch(circle)

        xl = self.remove_length_units(bounds[0])
        yl = self.remove_length_units(bounds[1])
        plt.xlim([-xl, 2*xl])
        plt.ylim([-yl, 2*yl])
        plt.draw()
        plt.pause(0.01)


# configs
DEFAULT_DIAMETER = 7.5 * DEFAULT_LENGTH_UNIT

def single_cell_config(config):
    # cell dimensions
    diameter = DEFAULT_DIAMETER
    volume = sphere_volume_from_diameter(diameter)
    velocity = config.get('velocity', 1 * DEFAULT_VELOCITY_UNIT)
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    location = config.get('location')
    if location:
        location = [loc * bounds[n] for n, loc in enumerate(location)]
    else:
        location = make_random_position(bounds)
    return {
        'boundary': {
            'location': location,
            'velocity': velocity,
            'volume': volume,
            'diameter': diameter,
            'mass': 5 * units.ng,
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


default_gd_config = {'bounds': DEFAULT_BOUNDS}
default_gd_config.update(cell_body_config({
    'bounds': DEFAULT_BOUNDS,
    'cell_ids': ['1', '2']}))

def test_growth_division(config=default_gd_config, settings={}):
    initial_cells_state = config['cells']

    # make the process
    multibody = Neighbors(config)
    experiment = process_in_experiment(multibody, settings)
    experiment.state._update_subschema(
        ('cells',), {
            'boundary': {
                'mass': {
                    '_updater': 'set',
                    '_divider': 'split'},
                'diameter': {
                    '_updater': 'set',
                    '_divider': 'split'},
                }})
    experiment.state._apply_subschemas()

    # get initial cell state
    experiment.state.set_value({'cells': initial_cells_state})
    cells_store = experiment.state.get_path(['cells'])

    ## run simulation
    # get simulation settings
    growth_rate = settings.get('growth_rate', 0.0006)
    growth_rate_noise = settings.get('growth_rate_noise', 0.0)
    division_volume = settings.get('division_volume', 0.4 * DEFAULT_LENGTH_UNIT ** 3)
    total_time = settings.get('total_time', 120)
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
                        'key': daughter_id,
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

            invoked_update.append((InvokeUpdate({'cells': update}), None))

        # update experiment
        experiment._send_updates(invoked_update)

    experiment.end()
    return experiment.emitter.get_data()

def multibody_neighbors_workflow(
        config={}, out_dir='out', filename='neighbors'):
    n_cells = 2
    cell_ids = [str(cell_id) for cell_id in range(n_cells)]

    bounds = DEFAULT_BOUNDS
    settings = {
        'growth_rate': 0.04,
        'growth_rate_noise': 0.02,
        'division_volume': sphere_volume_from_diameter(10 * DEFAULT_LENGTH_UNIT),
        'progress_bar': False,
        'display_info': False,
        'total_time': 100}
    gd_config = {
        'animate': True,
        'jitter_force': 1e-4,
        'bounds': bounds}
    body_config = {
        'velocity': 5 * DEFAULT_VELOCITY_UNIT,
        'bounds': bounds,
        'cell_ids': cell_ids}
    gd_config.update(cell_body_config(body_config))
    gd_data = test_growth_division(gd_config, settings)



if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    multibody_neighbors_workflow({'animate': True}, out_dir)
