import math
import random

from vivarium.library.units import units

# constants
PI = math.pi


def random_location(
        bounds,
        center=None,
        distance_from_center=None,
        excluded_distance_from_center=None
):
    """
    generate a single random location within `bounds`, and within `distance_from_center`
    of a provided `center`. `excluded_distance_from_center` is an additional parameter
    that leaves an empty region around the center point.
    """
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


DEFAULT_LENGTH_UNIT = units.um
DEFAULT_BOUNDS = [200 * DEFAULT_LENGTH_UNIT, 200 * DEFAULT_LENGTH_UNIT]
