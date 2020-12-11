import random
import math
import copy

from vivarium.library.units import units, Quantity

import pymunk


PI = math.pi



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

def random_direction(velocity):
    angle = random.uniform(0, 360)
    return (velocity * math.cos(angle), velocity * math.sin(angle))


class NullScreen(object):
    def update_screen(self):
        pass
    def configure(self, config):
        pass


class PymunkMinimal(object):
    """
    Multibody object for interfacing with pymunk
    """

    defaults = {
        'cell_shape': 'segment',
        # hardcoded parameters
        'elasticity': 0.9,
        'friction': 0.9,
        'physics_dt': 0.01,
        # configured parameters
        'length_unit': units.mm,
        'mass_unit': units.fg,
        'velocity_unit': units.mm / units.s,
        'bounds': [20, 20],
        'barriers': False,
        'initial_cells': {},
    }

    def __init__(self, config):
        # hardcoded parameters
        self.elasticity = self.defaults['elasticity']
        self.friction = self.defaults['friction']
        self.physics_dt = config.get('physics_dt', self.defaults['physics_dt'])

        # configured parameters
        self.length_unit = config.get('length_unit', self.defaults['length_unit'])
        self.mass_unit = config.get('mass_unit', self.defaults['mass_unit'])
        self.velocity_unit = config.get('velocity_unit', self.defaults['velocity_unit'])
        self.cell_shape = config.get('cell_shape', self.defaults['cell_shape'])
        bounds = config.get('bounds', self.defaults['bounds'])
        self.bounds = [bound.to(self.length_unit).magnitude for bound in bounds]
        barriers = config.get('barriers', self.defaults['barriers'])

        # initialize pymunk space
        self.space = pymunk.Space()

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
        self.space.step(timestep)


    def add_barriers(self, bounds, barriers):
        """ Create static barriers """
        thickness = 50.0
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
            self.space.add(line)

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

    def set_velocity(self, body_id, velocity):
        body, shape = self.bodies[body_id]
        if isinstance(velocity, Quantity):
            v = velocity.to(self.velocity_unit).magnitude
            body.velocity = random_direction(v)
        else:
            body.velocity = random_direction(velocity)

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
            pymunk_bodies[bodies_id]['boundary']['location'] = [loc.to(self.length_unit).magnitude for loc in specs['boundary']['location']]
            # convert diameter
            pymunk_bodies[bodies_id]['boundary']['diameter'] = specs['boundary']['diameter'].to(self.length_unit).magnitude
            # convert mass
            pymunk_bodies[bodies_id]['boundary']['mass'] = specs['boundary']['mass'].to(self.mass_unit).magnitude
        return pymunk_bodies


    def pymunk_location_to_body_units(self, bodies):
        for body_id, location in bodies.items():
            bodies[body_id] = [(loc * self.length_unit) for loc in location]
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



def test_minimal(
        total_time=2,
        cell_shape='rectangle',
        velocity=10.0 * units.um / units.min,  #  um/minute
        n_cells=1,
):

    bounds = [500 * units.um, 500 * units.um]
    center_location = [0.5*loc.magnitude for loc in bounds]
    cell_ids = [str(cell_idx) for cell_idx in range(n_cells)]
    cells = {
        cell_id: {
            'boundary': {
                'location': center_location,
                'volume': 15,
                'diameter': 30,
                'mass': 1,
                'velocity': velocity,
            }}
        for cell_id in cell_ids
    }
    config = {
        'cell_shape': cell_shape,
        'bounds': bounds,
        'barriers': False,
        'initial_cells': cells,
    }
    multibody = PymunkMinimal(config)

    # run simulation
    time = 0
    time_step = 1
    while time < total_time:
        time += time_step
        multibody.run(time_step)
        for cell_id in cell_ids:
            multibody.set_velocity(cell_id, velocity)

        print(multibody.get_body_positions())


if __name__ == '__main__':
    test_minimal(
        total_time=10)
