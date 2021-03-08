import random
import math

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
    angle = random.uniform(0, 2*PI)
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
        # configured parameters
        'physics_dt': 0.01,
        'bounds': [20, 20],
        'initial_cells': {},
    }

    def __init__(self, config):
        # hardcoded parameters
        self.elasticity = self.defaults['elasticity']
        self.friction = self.defaults['friction']

        # configured parameters
        self.physics_dt = config.get('physics_dt', self.defaults['physics_dt'])
        self.cell_shape = config.get('cell_shape', self.defaults['cell_shape'])
        bounds = config.get('bounds', self.defaults['bounds'])

        # initialize pymunk space
        self.space = pymunk.Space()

        # add static barriers
        self.add_barriers(self.space, bounds)

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


    def add_barriers(self, space, bounds):
        """ Create static barriers """
        static_body = space.static_body

        thickness = 500
        offset = thickness
        x0 = 0 - offset
        x1 = bounds[0] + offset
        y0 = 0 - offset
        y1 = bounds[1] + offset
        ps = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for i in range(4):
            segment = pymunk.Segment(static_body, ps[i], ps[(i+1) % 4], thickness)
            segment.elasticity = 0.0
            segment.friction = 0.8
            space.add(segment)

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
        # body.velocity = random_direction(velocity)

        # import ipdb; ipdb.set_trace()


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

        if 'velocity' in boundary:
            self.set_velocity(body_id, boundary['velocity'])

    def update_bodies(self, bodies):
        # if an cell has been removed from the cells store,
        # remove it from space and bodies

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
        return {
            body_id: self.get_body_position(body_id)
            for body_id in self.bodies.keys()}



def test_minimal(
        total_time=2,
        cell_shape='rectangle',
        velocity=10.0,  #  um / min
        n_cells=1,
):

    bounds = [500, 500]
    center_location = [0.5*loc for loc in bounds]
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
        total_time=100)
