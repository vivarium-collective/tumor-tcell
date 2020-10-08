from __future__ import absolute_import, division, print_function

import random
import math
import copy

from vivarium.library.units import units

# pymunk imports
import pymunkoptions
pymunkoptions.options["debug"] = False
import pymunk


PI = math.pi

DEBUG_SIZE = 600  # size of the pygame debug screen



def random_body_position(body):
    ''' pick a random point along the boundary'''
    diameter = body.diameter
    if random.randint(0, 1) == 0:
        # force along ends
        if random.randint(0, 1) == 0:
            # force on the left end
            location = (random.uniform(0, diameter), 0)
        else:
            # force on the right end
            location = (random.uniform(0, diameter), diameter)
    else:
        # force along length
        if random.randint(0, 1) == 0:
            # force on the bottom end
            location = (0, random.uniform(0, diameter))
        else:
            # force on the top end
            location = (diameter, random.uniform(0, diameter))
    return location


class NullScreen(object):
    def update_screen(self):
        pass
    def configure(self, config):
        pass


class PymunkMultibody(object):
    """
    Multibody object for interfacing with pymunk
    """

    defaults = {
        'cell_shape': 'segment',
        # hardcoded parameters
        'elasticity': 0.9,
        'damping': 0.5,  # 1 is no damping, 0 is full damping
        'angular_damping': 0.8,
        'friction': 0.9,  # does this do anything?
        'physics_dt': 0.001,
        'force_scaling': 1e2,  # scales from pN
        # configured parameters
        'pymunk_length_unit': units.mm,
        'pymunk_mass_unit': units.fg,
        'jitter_force': 1e-3,
        'bounds': [20, 20],
        'barriers': False,
        'initial_cells': {},
        # for debugging
        'screen': None,
    }

    def __init__(self, config):
        # hardcoded parameters
        self.elasticity = self.defaults['elasticity']
        self.friction = self.defaults['friction']
        self.damping = self.defaults['damping']
        self.angular_damping = self.defaults['angular_damping']
        self.physics_dt = config.get('physics_dt', self.defaults['physics_dt'])
        self.force_scaling = self.defaults['force_scaling']

        # configured parameters
        self.pymunk_length_unit = config.get('pymunk_length_unit', self.defaults['pymunk_length_unit'])
        self.pymunk_mass_unit = config.get('pymunk_mass_unit', self.defaults['pymunk_mass_unit'])
        self.cell_shape = config.get('cell_shape', self.defaults['cell_shape'])
        self.jitter_force = config.get('jitter_force', self.defaults['jitter_force'])
        bounds = config.get('bounds', self.defaults['bounds'])
        self.bounds = [bound.to(self.pymunk_length_unit).magnitude for bound in bounds]
        barriers = config.get('barriers', self.defaults['barriers'])

        # initialize pymunk space
        self.space = pymunk.Space()

        # debug screen
        self.screen = config.get('screen')
        if self.screen is None:
            self.screen = NullScreen()
        self.screen.configure({
            'space': self.space,
            'bounds': self.bounds})

        # add static barriers
        self.add_barriers(self.bounds, barriers)

        # initialize cells
        initial_cells = config.get('initial_cells', self.defaults['initial_cells'])
        self.bodies = {}
        for cell_id, specs in initial_cells.items():
            self.add_body_from_center(cell_id, specs)

    def run(self, timestep):
        if self.physics_dt > timestep:
            print('timestep skipped by pymunk_multibody: {}'.format(timestep))
            return

        time = 0
        while time < timestep:
            time += self.physics_dt

            # apply forces
            for body in self.space.bodies:
                self.apply_jitter_force(body)
                self.apply_viscous_force(body)

            # run for a physics timestep
            self.space.step(self.physics_dt)

        self.screen.update_screen()

    def apply_jitter_force(self, body):
        jitter_location = random_body_position(body)
        jitter_force = [
            random.normalvariate(0, self.jitter_force),
            random.normalvariate(0, self.jitter_force)]
        scaled_jitter_force = [
            force * self.force_scaling
            for force in jitter_force]
        body.apply_impulse_at_local_point(
            scaled_jitter_force,
            jitter_location)

    def apply_viscous_force(self, body):
        # dampen velocity
        body.velocity = body.velocity * self.damping + (body.force / body.mass) * self.physics_dt

    def add_barriers(self, bounds, barriers):
        """ Create static barriers """
        thickness = 5.0
        offset = thickness
        x_bound = bounds[0]
        y_bound = bounds[1]

        static_body = self.space.static_body
        static_lines = [
            pymunk.Segment(
                static_body,
                (0.0-offset, 0.0-offset),
                (x_bound+offset, 0.0-offset),
                thickness),
            pymunk.Segment(
                static_body,
                (x_bound+offset, 0.0-offset),
                (x_bound+offset, y_bound+offset),
                thickness),
            pymunk.Segment(
                static_body,
                (x_bound+offset, y_bound+offset),
                (0.0-offset, y_bound+offset),
                thickness),
            pymunk.Segment(
                static_body,
                (0.0-offset, y_bound+offset),
                (0.0-offset, 0.0-offset),
                thickness),
        ]

        if barriers:
            spacer_thickness = barriers.get('spacer_thickness', 0.1)
            channel_height = barriers.get('channel_height')
            channel_space = barriers.get('channel_space')
            n_lines = math.floor(x_bound/channel_space)

            machine_lines = [
                pymunk.Segment(
                    static_body,
                    (channel_space * line, 0),
                    (channel_space * line, channel_height), spacer_thickness)
                for line in range(n_lines)]
            static_lines += machine_lines

        for line in static_lines:
            line.elasticity = 0.0  # bounce
            line.friction = 0.8
        self.space.add(static_lines)

    def get_shape(self, boundary):
        '''
        shape documentation at: https://pymunk-tutorial.readthedocs.io/en/latest/shape/shape.html
        '''

        diameter = boundary['diameter']
        half_diameter = diameter / 2
        shape = pymunk.Circle(None, radius=half_diameter, offset=(0, 0))

        return shape

    def get_inertia(self, shape, mass):
        radius = shape.radius
        inertia = pymunk.moment_for_circle(mass, radius, radius)
        return inertia

    def add_body_from_center(self, body_id, specs):
        boundary = specs['boundary']
        mass = boundary['mass']
        center_position = boundary['location']
        diameter = boundary['diameter']

        # get shape, inertia, make body, assign body to shape
        shape = self.get_shape(boundary)
        inertia = self.get_inertia(shape, mass)
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = (
            center_position[0],
            center_position[1])
        body.diameter = diameter

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add body to cells dictionary
        self.bodies[body_id] = (body, shape)

    def update_body(self, body_id, specs):
        boundary = specs['boundary']
        diameter = boundary['diameter']
        mass = boundary['mass']

        body, shape = self.bodies[body_id]
        position = body.position

        # get shape, inertia, make body, assign body to shape
        new_shape = self.get_shape(boundary)
        inertia = self.get_inertia(new_shape, mass)
        new_body = pymunk.Body(mass, inertia)
        new_shape.body = new_body

        new_body.position = position
        new_body.velocity = body.velocity
        new_body.angular_velocity = body.angular_velocity
        new_body.diameter = diameter

        new_shape.elasticity = shape.elasticity
        new_shape.friction = shape.friction

        # swap bodies
        self.space.remove(body, shape)
        self.space.add(new_body, new_shape)

        # update body
        self.bodies[body_id] = (new_body, new_shape)


    def bodies_to_pymunk_units(self, bodies):
        pymunk_bodies = copy.deepcopy(bodies)
        for bodies_id, specs in pymunk_bodies.items():
            # convert location
            pymunk_bodies[bodies_id]['boundary']['location'] = [loc.to(self.pymunk_length_unit).magnitude for loc in specs['boundary']['location']]
            # convert diameter
            pymunk_bodies[bodies_id]['boundary']['diameter'] = specs['boundary']['diameter'].to(self.pymunk_length_unit).magnitude
            # convert mass
            pymunk_bodies[bodies_id]['boundary']['mass'] = specs['boundary']['mass'].to(self.pymunk_mass_unit).magnitude
        return pymunk_bodies


    def pymunk_location_to_body_units(self, bodies):
        for body_id, location in bodies.items():
            bodies[body_id] = [(loc * self.pymunk_length_unit) for loc in location]
        return bodies


    def update_bodies(self, raw_bodies):
        # if an cell has been removed from the cells store,
        # remove it from space and bodies

        # convert to pymunk_units
        bodies = self.bodies_to_pymunk_units(raw_bodies)

        removed_bodies = [
            body_id for body_id in self.bodies.keys()
            if body_id not in bodies.keys()]
        for body_id in removed_bodies:
            body, shape = self.bodies[body_id]
            self.space.remove(body, shape)
            del self.bodies[body_id]

        # update cells, add new cells
        for body_id, specs in bodies.items():
            if body_id in self.bodies:
                self.update_body(body_id, specs)
            else:
                self.add_body_from_center(body_id, specs)

    def get_body_position(self, cell_id):
        body, shape = self.bodies[cell_id]
        return tuple(pos for pos in body.position)

    def get_body_positions(self):
        bodies = {
            body_id: self.get_body_position(body_id)
            for body_id in self.bodies.keys()}
        return self.pymunk_location_to_body_units(bodies)



def test_multibody(
        total_time=2,
        cell_shape='rectangle',
        n_cells=1,
        jitter_force=1e1,
        screen=None):

    bounds = [500 * units.um, 500 * units.um]
    center_location = [0.5*loc.magnitude for loc in bounds]
    cells = {
        str(cell_idx): {
            'boundary': {
                'location': center_location,
                'volume': 15,
                'diameter': 30,
                'mass': 1}}
        for cell_idx in range(n_cells)
    }
    config = {
        'cell_shape': cell_shape,
        'jitter_force': jitter_force,
        'bounds': bounds,
        'barriers': False,
        'initial_cells': cells,
        'screen': screen
    }
    multibody = PymunkMultibody(config)

    # run simulation
    time = 0
    time_step = 1
    while time < total_time:
        time += time_step
        multibody.run(time_step)


if __name__ == '__main__':
    test_multibody(10)
