"""
Crystallographic viewing directions
Triclinic: -
Monoclinic: b(c)
Orthorhombic: a b c
Tetragonal: c a [110]
Trigonal: c a [210]
Hexagonal: c a [210]
cubic: a [111] [110]
"""
from functions import *
from readinginput import *

def sym_to_lattice(sym_no):
    if 1 <= sym_no <= 2:
        lattice = 'Triclinic'
    elif 3 <= sym_no <= 15:
        lattice = 'Monoclinic'
    elif 16 <= sym_no <= 74:
        lattice = 'Orthorhombic'
    elif 75 <= sym_no <= 142:
        lattice = 'Tetragonal'
    elif 143 <= sym_no <= 167:
        lattice = 'Trigonal'
    elif 168 <= sym_no <= 194:
        lattice = 'Hexagonal'
    elif 195 <= sym_no <= 230:
        lattice = 'Cubic'
    return lattice


# initial wyckoff_positions are direct coordinates
# initial_positions are car coordinaties
def move_rigid_to_wyckoff(initial_positions, wyckoff_positions, rigid_type, cell):
    wyckoff_move = np.array([0.0, 0.0, 0.0])
    random_wyckoff = np.array([0.0, 0.0, 0.0])
    rigid_class = rigid_select(rigid_type)(1, 1.0)
    if wyckoff_positions[0] != 'x':
        random_wyckoff[0] = wyckoff_positions[0]
    elif wyckoff_positions[0] == 'x':
        random_wyckoff[0] = round(random.uniform(0, 1), 3)
    if wyckoff_positions[1] != 'y':
        random_wyckoff[1] = wyckoff_positions[1]
    elif wyckoff_positions[1] == 'y':
        random_wyckoff[1] = round(random.uniform(0, 1), 3)
    if wyckoff_positions[2] != 'z':
        random_wyckoff[2] = wyckoff_positions[2]
    elif wyckoff_positions[2] == 'z':
        random_wyckoff[2] = round(random.uniform(0, 1), 3)
    random_wyckoff_car = direct_cartesian_transform(random_wyckoff, cell, 'DtoC')
    rigid_center = rigid_class.center(initial_positions)
    wyckoff_move = random_wyckoff_car - rigid_center
    final_positions = initial_positions + wyckoff_move
    return final_positions
    

def rigid_select(rigid_type):
    if rigid_type == 'Tetrahedron':
        rigid_sele = Tetrahedron
    if rigid_type == 'GeSe':
        rigid_sele = GeSe
    return rigid_sele


class rigid_rotation:
    def __init__(self, positions, center):
        self.positions = positions
        self.center = center

    def random_rotation(self, sym_no, sym_element, angle_range, cell):
        lattice_type = sym_to_lattice(sym_no)
        random_angle = round(random.uniform(-angle_range, angle_range), 0)
        random_vector = np.random.rand(3)
        random_unit_vector = random_vector / np.linalg.norm(random_vector)
        r_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if lattice_type == 'Triclinic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
        elif lattice_type == 'Monoclinic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
            elif sym_element == ['2'] or ['m'] or ['2/m']:
                b_unit = cell[1] / np.linalg.norm(cell[1])
                r_matrix = rotation_matrix(b_unit, random_angle)
        elif lattice_type == 'Orthorhombic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
            elif sym_element == ['x', 'y', '2'] or sym_element == ['x', 'y', 'm'] or sym_element == ['x', 'y', '2/m']:
                r_matrix = rotation_matrix([0, 0, 1], random_angle)
            elif sym_element == ['x', '2', 'z'] or sym_element == ['x', 'm', 'z'] or sym_element == ['x', '2/m', 'z']:
                r_matrix = rotation_matrix([0, 1, 0], random_angle)
            elif sym_element == ['2', 'y', 'z'] or sym_element == ['m', 'y', 'z'] or sym_element == ['2/m', 'y', 'z']:
                r_matrix = rotation_matrix([1, 0, 0], random_angle)
            elif sym_element == ['2', '2', '2']:
                ra = rotation_matrix([1, 0, 0], random_angle)
                rb = rotation_matrix([0, 1, 0], random_angle)
                rc = rotation_matrix([0, 0, 1], random_angle)
                r_matrix = ra @ rb @ rc
       
        move_array = np.array([0, 0, 0]) - self.center
        origin_positions = self.positions + move_array
        rotate_positions = origin_positions @ r_matrix
        final_positions = rotate_positions - move_array
        return final_positions


def setup_random_rotation(sym_no, sym_element, angle_range, cell):
        lattice_type = sym_to_lattice(sym_no)
        random_angle = round(random.uniform(-angle_range, angle_range), 0)
        random_vector = np.random.rand(3)
        random_unit_vector = random_vector / np.linalg.norm(random_vector)
        # r_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if lattice_type == 'Triclinic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
        elif lattice_type == 'Monoclinic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
            elif sym_element == ['2'] or ['m'] or ['2/m']:
                b_unit = cell[1] / np.linalg.norm(cell[1])
                r_matrix = rotation_matrix(b_unit, random_angle)
        elif lattice_type == 'Orthorhombic':
            if sym_element == ['1'] or sym_element == ['-1']:
                r_matrix = rotation_matrix(random_unit_vector, random_angle)
            elif sym_element == ['x', 'y', '2'] or sym_element == ['x', 'y', 'm'] or sym_element == ['x', 'y', '2/m']:
                r_matrix = rotation_matrix([0, 0, 1], random_angle)
            elif sym_element == ['x', '2', 'z'] or sym_element == ['x', 'm', 'z'] or sym_element == ['x', '2/m', 'z']:
                r_matrix = rotation_matrix([0, 1, 0], random_angle)
            elif sym_element == ['2', 'y', 'z'] or sym_element == ['m', 'y', 'z'] or sym_element == ['2/m', 'y', 'z']:
                r_matrix = rotation_matrix([1, 0, 0], random_angle)
            elif sym_element == ['2', '2', '2']:
                ra = rotation_matrix([1, 0, 0], random_angle)
                rb = rotation_matrix([0, 1, 0], random_angle)
                rc = rotation_matrix([0, 0, 1], random_angle)
                r_matrix = ra @ rb @ rc
        return r_matrix


class Tetrahedron:
    def __init__(self, sym_no, bond_len):
        self.sym_no = sym_no
        self.bond_len = bond_len
        self.atom_no = 5

    def atoms_position(self):
        site_ori = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]])
        site_fin = (site_ori / 1.73205) * self.bond_len
        return site_fin
    
    def center(self, atom_positions):
        center_position = atom_positions[0]
        return center_position

    def wyckoff_sym(self):
        lattice_type = sym_to_lattice(self.sym_no)
        if lattice_type == 'Triclinic':
            wyck_sym = [['1'], ['m'], ['2']]
        elif lattice_type == 'Monoclinic':
            wyck_sym = [['1'], ['m'], ['2']]
        elif lattice_type == 'Orthorhombic':
            wyck_sym = [['1'], ['m', 'y', 'z'], ['x', 'm', 'z'], ['x', 'y', 'm'], ['2', 'y', 'z'], 
                        ['x', '2', 'z'],['x', 'y', '2'], ['2', 'm', 'm'], ['m', '2', 'm'], ['m', 'm', '2']]
        elif lattice_type == 'Tetragonal':
            wyck_sym = [['1'], ['m']]
        elif lattice_type == 'Trigonal':
            wyck_sym = [['1'], ['m']]
        elif lattice_type == 'Hexagonal':
            wyck_sym = [['1'], ['2', 'y', 'z'], ['x', '2', 'z'], ['x', 'y', '2']]
        return wyck_sym


    def rotaiton_matrix(self, sym_element):
        lattice_type = sym_to_lattice(self.sym_no)
        r_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if lattice_type == 'Monoclinic':
            if sym_element == ['m']:
                r_matrix = rotation_matrix([1, 0, 0], 45)
            elif sym_element == ['2']:
                r_matrix = rotation_matrix([1, 0, 0], 45) @ rotation_matrix([0, 0, 1], 90)
        elif lattice_type == 'Orthorhombic':
            if sym_element == ['m', 'y', 'z'] or sym_element == ['x', 'm', 'z'] or \
               sym_element == ['x', 'y', '2'] or sym_element == ['m', 'm', '2']:
                r_matrix = rotation_matrix([0, 0, 1], 45)
            elif sym_element == ['x', 'y', 'm'] or sym_element == ['x', '2', 'z'] or \
                 sym_element == ['m', '2', 'm']:
                r_matrix == rotation_matrix([0, 0, 1], 45) @ rotation_matrix([1, 0, 0], 90)
            elif sym_element == ['2', 'y', 'z'] or sym_element == ['2', 'm', 'm']:
                r_matrix == rotation_matrix([0, 0, 1], 45) @ rotation_matrix([0, 1, 0], 90)
        elif lattice_type == 'Tetragonal':
            pass
        elif lattice_type == 'Trigonal':
            pass    
        elif lattice_type == 'Hexagonal':
            pass
        elif lattice_type == 'cubic':
            pass         
        return r_matrix

class GeSe:
    def __init__(self, sym_no, bond_len):
        self.sym_no = sym_no
        self.bond_len = bond_len
        self.atom_no = 8

    def atoms_position(self):
        site_ori = np.array([[0.00000000e+00, 0.00000000e+00, 8.66025000e-01],
                             [0.00000000e+00, 0.00000000e+00, -8.66025000e-01],
                             [1.41421356e+00, 8.16496477e-01, 1.44337542e+00],
                             [1.11022302e-16, -1.63299327e+00, 1.44337498e+00],
                             [-1.41421356e+00, 8.16496477e-01, 1.44337542e+00],
                             [-8.88175036e-08, 1.63299311e+00, -1.44337542e+00],
                             [-1.41421366e+00, -8.16496635e-01, -1.44337498e+00],
                             [1.41421347e+00, -8.16496631e-01, -1.44337542e+00]])
        site_fin = (site_ori / 1.73205) * self.bond_len
        return site_fin

    def center(self, atom_positions):
        center_position = (atom_positions[0] + atom_positions[1]) / 2.0
        return center_position

    def wyckoff_sym(self):
        lattice_type = sym_to_lattice(self.sym_no)
        if lattice_type == 'Triclinic':
            wyck_sym = [['1'], ['-1']]
        elif lattice_type == 'Monoclinic':
            wyck_sym = [['1'], ['-1'], ['m']]
        elif lattice_type == 'Orthorhombic':
            wyck_sym = [['1'], ['m', 'y', 'z'], ['x', 'm', 'z'], ['x', 'y', 'm'],]
        elif lattice_type == 'Tetragonal':
            wyck_sym = [['1'], ['m']]
        elif lattice_type == 'Trigonal':
            wyck_sym = [['1'], ['m']]
        elif lattice_type == 'Hexagonal':
            wyck_sym = [['1']]
        return wyck_sym

    def rotaiton_matrix(self, sym_element):
        lattice_type = sym_to_lattice(self.sym_no)
        r_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if lattice_type == 'Monoclinic':
            if sym_element == ['m']:
                r_matrix = rotation_matrix([0, 0, 1], 90)
        elif lattice_type == 'Orthorhombic':
            if sym_element == ['x', 'm', 'z']:
                r_matrix = rotation_matrix([0, 0, 1], 90)
            elif sym_element == ['x', 'y', 'm']:
                r_matrix == rotation_matrix([1, 0, 0], 90) @ rotation_matrix([0, 1, 0], 90)
        elif lattice_type == 'Tetragonal':
            pass
        elif lattice_type == 'Trigonal':
            pass    
        elif lattice_type == 'Hexagonal':
            pass
        elif lattice_type == 'cubic':
            pass         
        return r_matrix
