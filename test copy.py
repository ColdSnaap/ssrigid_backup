import numpy as np
from functions import direct_cartesian_transform, move_atoms_into_box

cell = np.array([[6.9336127435, 0., 0.],
        [0., 6.9336127435, 0.],
        [0., 0., 6.9336127435]])

position = np.array([3., 3., 3.])
add_matrix = np.array([0., 0.5, 0.])
# direct coor plus

# position_dir = direct_cartesian_transform(position, cell, 'CtoD')
# print(position_dir)
def dir_coor_add(cell, position, add_matrix):
    