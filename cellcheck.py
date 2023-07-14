# from readinginput import *
import numpy as np
import math
from functions import distance

def angle_of_vectors(vector1, vector2):
    dotproduct = np.dot(vector1, vector2)
    modofvector1 = math.sqrt(vector1[0]**2 + vector1[1]**2 + vector1[2]**2)
    modofvector2 = math.sqrt(vector2[0]**2 + vector2[1]**2 + vector2[2]**2)
    angle = dotproduct / (modofvector1 * modofvector2)
    angleindegree = math.degrees(math.acos(angle))
    return angleindegree


def lattice_angle(lattice):
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    normal_a_b = np.cross(a, b)
    angle = angle_of_vectors(normal_a_b, c)
    return angle


def rigid_check(rigid_type, positions):
    pass