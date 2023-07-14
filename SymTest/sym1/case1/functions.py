import math
import os
import random
import shutil
from tkinter import Image
import numpy as np
from collections import Counter
from ase.spacegroup import crystal


def step_algorithm(algorithm, temp):
    shift_current = 0
    if algorithm == 'normal':
        s = np.random.normal(0, temp, 1)
        shift_current_1 = float(s[0])
        shift_current = 0.3 * shift_current_1
    return shift_current


def temp_algorithm(algorithm, temp_ini, script_count):
    current_temp = 0
    if algorithm == 'normal':
        current_temp = temp_ini / math.log(script_count)
    if algorithm == 'fast':
        # current_temp = temp_ini * (0.99 ** script_count)
        current_temp = 1.5
    return current_temp


def initial_cell(cell_type, atom_amount, density):
    if cell_type == 'cube':
            cube_len = round((atom_amount / density) ** (1./3), 10)
            cell = np.zeros((3, 3))
            for i in range(3):
                cell[i][i] = cube_len
    elif cell_type == 'hex':
        hex_len = round(((atom_amount / density) ** (1 / 3.0)) / math.sqrt(3.0), 10)
        cell = np.zeros((3, 3))
        cell[0][0], cell[2][2] = hex_len, hex_len
        cell[1][0] = - (hex_len / 2.0)
        cell[1][1] = (math.sqrt(3.0) / 2.0) * hex_len
    return cell


def direct_cartesian_transform(atom_positions, cell, transform_direction):
    atom_positions_final = atom_positions.copy()
    a1, a2, a3 = cell[0], cell[1], cell[2]
    # Build the matrix of lattice vectors stored column-wise
    # and get its inverse
    a_array = np.vstack([a1, a2, a3]).T
    a_array_inv = np.linalg.inv(a_array)
    if transform_direction == 'CtoD':
        atom_positions_final = np.matmul(a_array_inv, atom_positions.T).T
    elif transform_direction == 'DtoC':
        atom_positions_final = np.matmul(a_array, atom_positions.T).T
    return atom_positions_final


# atom positions are car coordinates
def move_atoms_into_box(atom_positions, cell):
    atom_amount = len(atom_positions)
    atom_positions_dir = direct_cartesian_transform(atom_positions, cell, 'CtoD')
    for i in range(atom_amount):
        for j in range(3):
            if atom_positions_dir[i][j] >= 1.0:
                while atom_positions_dir[i][j] >= 1.0:
                    atom_positions_dir[i][j] -= 1.0
            elif atom_positions_dir[i][j] < 0.0:
                while atom_positions_dir[i][j] < 0.0:
                    atom_positions_dir[i][j] += 1.0
    atom_positions_car = direct_cartesian_transform(atom_positions_dir, cell, 'DtoC')
    return atom_positions_car


def rotation_matrix(axis, angle):
    # make sure axis is unit vector
    if axis[0]**2 + axis[1]**2 + axis[2]**2 != 1.0:
        x = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        axis_after_x = round(axis[0] / x, 6)
        axis_after_y = round(axis[1] / x, 6)
        axis_after_z = round(axis[2] / x, 6)
        axis_unit = [axis_after_x, axis_after_y, axis_after_z]
    else:
        axis_unit = axis
    
    y = float(angle / 360) * (2.0 * math.pi)
    r11 = math.cos(y) + (axis_unit[0]**2) * (1 - math.cos(y))
    r12 = axis_unit[0] * axis_unit[1] * (1 - math.cos(y)) - (axis_unit[2] * math.sin(y))
    r13 = axis_unit[0] * axis_unit[2] * (1 - math.cos(y)) + axis_unit[1] * math.sin(y)
    r21 = axis_unit[1] * axis_unit[0] * (1 - math.cos(y)) + axis_unit[2] * math.sin(y)
    r22 = math.cos(y) + axis_unit[1]**2 * (1 - math.cos(y))
    r23 = axis_unit[1] * axis_unit[2] * (1 - math.cos(y)) - axis_unit[0] * math.sin(y)
    r31 = axis_unit[2] * axis_unit[0] * (1 - math.cos(y)) - axis_unit[1] * math.sin(y)
    r32 = axis_unit[2] * axis_unit[1] * (1 - math.cos(y)) + axis_unit[0] * math.sin(y)
    r33 = math.cos(y) + axis_unit[2]**2 * (1 - math.cos(y))

    r_matrix = np.array([[r11, r12, r13],
                         [r21, r22, r23],
                         [r31, r32, r33]])
    return r_matrix


def sym_information(sym_no, sym_dir, case_no):
    case_dir = sym_dir[case_no]
    rigid_list = case_dir['rigid']
    single_list = case_dir['single']
    with open('sym_information', 'w') as sym_info:
        sym_info.write('SpaceGroup : ' + f'{sym_no}' + '\n')
        sym_info.write('WyckoffRigid : ')
        for i in rigid_list:
            sym_info.write(f'{i} ')
        sym_info.write('\n')
        sym_info.write('WyckoffSingle : ')
        for i in single_list:
            sym_info.write(f'{i} ')


class RigidSingleNo:
    def __init__(self):
        with open('sym_information', 'r') as test:
            line = test.readlines()        
            # get rigid_no and single_no        
            count = 0
            for i in line:
                if 'WyckoffRigid' in i:
                    line_split = line[count].split()
                    # print(line_split)
                    for j in range(len(line_split)):
                        if line_split[j] == ':':
                            list_final1 = line_split[j+1:]
                elif 'WyckoffSingle' in i:
                    line_split = line[count].split()
                    # print(line_split)
                    for j in range(len(line_split)):
                        if line_split[j] == ':':
                            list_final2 = line_split[j+1:]
                count += 1
            list_final_rigid = list_final1.copy()
            list_final_single = list_final2.copy()
            self.rigid_no = len(list_final_rigid)
            self.single_no = len(list_final_single)


def contcar_to_rigid(contcar_position):
    get_no = RigidSingleNo()
    rigid_no = get_no.rigid_no
    single_no = get_no.single_no
    with open('sym_information', 'r') as test:
        line = test.readlines()        
        for n in range(1, rigid_no + 1):
            rigid_array = np.array([[.1, .1, .1]])
            count = 0
            for i in line:
                if f'rigid_{n}' in i:
                    # print(f'yes, the line is {count}')
                    line_split = line[count].split()
                    # print(line_split)
                    for j in range(len(line_split)):
                        if line_split[j] == ':':
                            list_final1 = line_split[j+1:]
                count += 1
            list_final = list_final1.copy()
            for i in range(len(list_final1)):
                list_final[i] = int(list_final1[i])
            
            for i in list_final:
                rigid_array = np.append(rigid_array, np.array([contcar_position[i]]), axis=0)
            
            rigid_array_f = rigid_array[1:]
            with open(f'rigid_{n}', 'r') as rigid_read:
                rigid_read_line = rigid_read.readlines()
            with open(f'rigid_{n}', 'w') as rigid_write:
                rigid_write.write(rigid_read_line[0])
                for j in range(len(rigid_array_f)):
                    c_x = rigid_array_f[j][0]
                    c_y = rigid_array_f[j][1]
                    c_z = rigid_array_f[j][2]
                    rigid_write.write(f'{c_x} {c_y} {c_z}' + '\n')
        
        for n in range(1, single_no + 1):
            rigid_array = np.array([[.1, .1, .1]])
            count = 0
            for i in line:
                if f'single_{n}' in i:
                    # print(f'yes, the line is {count}')
                    line_split = line[count].split()
                    # print(line_split)
                    for j in range(len(line_split)):
                        if line_split[j] == ':':
                            list_final1 = line_split[j+1:]
                    break
                count += 1
            list_final = list_final1.copy()
            for i in range(len(list_final1)):
                list_final[i] = int(list_final1[i])
            
            for i in list_final:
                rigid_array = np.append(rigid_array, np.array([contcar_position[i]]), axis=0)
            
            rigid_array_f = rigid_array[1:]
            with open(f'single_{n}', 'r') as rigid_read:
                rigid_read_line = rigid_read.readlines()
            with open(f'single_{n}', 'w') as rigid_write:
                rigid_write.write(rigid_read_line[0])
                for j in range(len(rigid_array_f)):
                    c_x = rigid_array_f[j][0]
                    c_y = rigid_array_f[j][1]
                    c_z = rigid_array_f[j][2]
                    rigid_write.write(f'{c_x} {c_y} {c_z}' + '\n')


def generate_single_atoms(wyckoff_positions, cell):
    atom_positions = np.array([0.0, 0.0, 0.0])
    if wyckoff_positions[0] != 'x':
        atom_positions[0] = wyckoff_positions[0]
    else:
        atom_positions[0] = round(random.uniform(0, 1), 3)
    
    if wyckoff_positions[1] != 'y':
        atom_positions[1] = wyckoff_positions[1]
    else:
        atom_positions[1] = round(random.uniform(0, 1), 3)

    if wyckoff_positions[2] != 'z':
        atom_positions[2] = wyckoff_positions[2]
    else:
        atom_positions[2] = round(random.uniform(0, 1), 3)

    atom_positions = direct_cartesian_transform(atom_positions, cell, 'DtoC')
    return atom_positions


# from ['1', '2', '3'] to np.array([1, 2, 3])
def number_in_list_trans(string):
    # print(f'string:{string}')
    output = np.array([.0, .0, .0])
    for i in range(3):
        x = float(string[i])
        output[i] = x
    # print(f'output:{output}')
    return output


def same_atom_check(atom1, atom2, cell):
    check = 0
    atom1_image = AtomImage(atom1, cell)
    close_dis = atom1_image.close_image_dis(atom2)
    if close_dis < 0.001:
        check = 1
    return check


def write_rigid_trans(file_name, sg_no, cell):
    atom_list_final = []
    atom_cor_final = np.array([[.1, .1, .1]])
    with open(file_name, 'r') as rigid:
        line = rigid.readlines()
        rigid_list = line[0].split()
        for j in range(len(rigid_list)):
            atom_name = rigid_list[j]
            atom_list = [atom_name]
            list_j = line[j+1].split()
            atom_array = number_in_list_trans(list_j)
            atom_array_dir = direct_cartesian_transform(atom_array, cell, 'CtoD')
            atom_cor = [(atom_array_dir[0], atom_array_dir[1], atom_array_dir[2])]
            trans1 = crystal(atom_name, atom_cor, spacegroup=sg_no, cell=cell)
            trans = trans1.positions
            len_trans = len(trans)
            atom_list_trans = atom_list * len_trans
            atom_list_final.extend(atom_list_trans)
            atom_cor_final = np.append(atom_cor_final, trans, axis=0)
    the_list = atom_list_final
    the_coor = atom_cor_final[1:]
    # remove same atom
    remove_list = []
    for i in range(len(the_coor)-1):
        if i not in remove_list:
            for j in range(i+1, len(the_coor)):
                check_same_atom = same_atom_check(the_coor[i], the_coor[j], cell)
                if check_same_atom == 1:
                    remove_list.append(j)
    remove_list.sort()
    count = 0
    for i in remove_list:
        the_list.remove(the_list[i-count])
        the_coor = np.delete(the_coor, i-count, 0)
        count += 1

    with open(f'{file_name}_trans','w') as trans_write:
        for i in the_list:
            trans_write.write(f'{i} ')
        trans_write.write('\n')
        for i in range(len(the_coor)):
            trans_write.write(f'{the_coor[i][0]} {the_coor[i][1]} {the_coor[i][2]}' + '\n')


# from rigid_1 ... single_1 file to write a dictory consist of all elements and
# their positions
def write_poscar_dic():
    get_no = RigidSingleNo()
    rigid_no = get_no.rigid_no
    single_no = get_no.single_no
    dic = {} 
    for i in range(1, rigid_no + 1):
        with open(f'rigid_{i}_trans', 'r') as rigid:
            line = rigid.readlines()
            rigid_list = line[0].split()
            for j in range(len(rigid_list)):
                list_j = line[j+1].split()
                list_j_float = number_in_list_trans(list_j)
                if rigid_list[j] in dic:
                    dic[rigid_list[j]] = np.append(dic[rigid_list[j]], np.array([list_j_float]), axis=0)
                else:
                    dic[rigid_list[j]] = np.array([list_j_float])
    for i in range(1, single_no + 1):
        with open(f'single_{i}_trans', 'r') as single:
            line2 = single.readlines()
            single_list = line2[0].split()
            for j in range(len(single_list)):
                list_j = line2[j+1].split()
                list_j_float = number_in_list_trans(list_j)
                if single_list[j] in dic:
                    dic[single_list[j]] = np.append(dic[single_list[j]], np.array([list_j_float]), axis=0)
                else:
                    dic[single_list[j]] = np.array([list_j_float])
    return dic


# from the dictory you got from write_poscar_dic to write a poscar out of that
def dir_to_poscar(dir, sim_box):
    with open('POSCAR', 'w') as poscar:
        poscar.write('sym1.0\n')
        poscar.write('1.0\n')
        for i in range(3):
            poscar.write(str(sim_box[i][0]) + '  ' + str(sim_box[i][1]) + '  ' + str(sim_box[i][2]) + '\n')
        dir_key = dir.keys()
        for i in dir_key:
            poscar.write(f'{i} ')
        poscar.write('\n')
        for i in dir_key:
            ele_no = len(dir[i])
            poscar.write(f'{str(ele_no)} ')
        poscar.write('\n')
        poscar.write('Car\n')
        positions_ini = np.array([[.1, .1, .1]])
        for i in dir_key:
            positions_ini = np.concatenate((positions_ini, dir[i]), axis=0)
        positions = positions_ini[1:]
        positions_len = len(positions)
        for i in range(positions_len):
            x = round(positions[i][0], 8)
            y = round(positions[i][1], 8)
            z = round(positions[i][2], 8)
            poscar.write(str(x) + '  ' + str(y) + '  ' + str(z) + '\n')


# wirte atom index number in poscar
def write_atom_index(file, cell):
    index = []
    poscar_read = ReadPOSCAR('POSCAR')
    poscar_position = poscar_read.positions
    with open(f'{file}', 'r') as rigid:
        line = rigid.readlines()
        rigid_list = line[0].split()
        for i in range(len(rigid_list)):
            list_i = line[i+1].split()
            atom_array = number_in_list_trans(list_i)
            for j in range(len(poscar_position)):
                atom_image = AtomImage(poscar_position[j], cell)
                atom_dis = atom_image.close_image_dis(atom_array)
                if atom_dis < 0.0001:
                    index.append(j)
    with open('sym_information', 'a') as sym_info:
        sym_info.write('\n')
        sym_info.write(f'{file} : ')
        for i in index:
            sym_info.write(f'{i} ')


# from rigid_1 ... single_1 ... files to write a array for mc
def rigid_to_array(file_name):
    array = np.array([[.1, .1, .1]])
    with open(file_name, 'r') as file:
        line = file.readlines()
        rigid_list = line[0].split()
        for j in range(len(rigid_list)):
            list_j = line[j+1].split()
            atom_array = number_in_list_trans(list_j)
            array = np.append(array, [atom_array], axis=0)
        array_f = array[1:]
    return array_f


def array_to_rigid(array, file_name):
    with open(file_name, 'r') as file_read:
        line = file_read.readlines()
    with open(file_name, 'w') as file_write:
        file_write.write(line[0])
        for j in range(len(array)):
            c_x = array[j][0]
            c_y = array[j][1]
            c_z = array[j][2]
            file_write.write(f'{c_x} {c_y} {c_z}' + '\n')


class ReadPOSCAR:
    def __init__(self, file):
        with open(file, 'r') as ps:
            line = ps.readlines()
            lattice = np.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    lattice[i][j] = float(line[i + 2].split()[j])
            atom = line[6].split()
            atom_type = len(atom)
            atom_number = 0
            for i in range(atom_type):
                atom_number = atom_number + int(atom[i])
            atom_coordinate = np.zeros((atom_number, 3))
            for i in range(atom_number):
                for j in range(3):
                    atom_coordinate[i][j] = float(line[i + 8].split()[j])
        self.lattice = lattice
        self.atom_number = atom_number
        self.positions = atom_coordinate


def read_system_initial(directory, total_atom):
    os.chdir(directory)
    with open('system_initial', 'r') as sys:
        line = sys.readlines()
        atom_coordinate = np.zeros((total_atom, 3))
        for i in range(total_atom):
            for j in range(3):
                atom_coordinate[i][j] = float(line[i].split()[j])
    return atom_coordinate


def distance(atom1, atom2):
    x_dis = atom1[0] - atom2[0]
    y_dis = atom1[1] - atom2[1]
    z_dis = atom1[2] - atom2[2]
    dis = math.sqrt(x_dis ** 2 + y_dis ** 2 + z_dis ** 2)
    return dis


# atom_with_image is car coordinates
class AtomImage:
    def __init__(self, atom_with_image, cell):
        self.atom_with_image = atom_with_image
        self.cell = cell

    def image(self):    
        atom_with_image_dir = direct_cartesian_transform(self.atom_with_image, self.cell, 'CtoD')
        image = np.zeros((27, 3))
        for i in range(27):
            image[i] = atom_with_image_dir
        # images
        # change only 1 cord
        image[1][0] = atom_with_image_dir[0] + 1.0
        image[2][0] = atom_with_image_dir[0] - 1.0
        image[3][1] = atom_with_image_dir[1] + 1.0
        image[4][1] = atom_with_image_dir[1] - 1.0
        image[5][2] = atom_with_image_dir[2] + 1.0
        image[6][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (+x, +y) (+x, +z) (+y, +z)
        image[7][0] = atom_with_image_dir[0] + 1.0
        image[7][1] = atom_with_image_dir[1] + 1.0
        image[8][0] = atom_with_image_dir[0] + 1.0
        image[8][2] = atom_with_image_dir[2] + 1.0
        image[9][1] = atom_with_image_dir[1] + 1.0
        image[9][2] = atom_with_image_dir[2] + 1.0

        # change 2 cord (-x, -y) (-x, -z) (-y, -z)
        image[10][0] = atom_with_image_dir[0] - 1.0
        image[10][1] = atom_with_image_dir[1] - 1.0
        image[11][0] = atom_with_image_dir[0] - 1.0
        image[11][2] = atom_with_image_dir[2] - 1.0
        image[12][1] = atom_with_image_dir[1] - 1.0
        image[12][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (+x, -y) (+x, -z) (+y, -z)
        image[13][0] = atom_with_image_dir[0] + 1.0
        image[13][1] = atom_with_image_dir[1] - 1.0
        image[14][0] = atom_with_image_dir[0] + 1.0
        image[14][2] = atom_with_image_dir[2] - 1.0
        image[15][1] = atom_with_image_dir[1] + 1.0
        image[15][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (-x, +y) (-x, +z) (-y, +z)
        image[16][0] = atom_with_image_dir[0] - 1.0
        image[16][1] = atom_with_image_dir[1] + 1.0
        image[17][0] = atom_with_image_dir[0] - 1.0
        image[17][2] = atom_with_image_dir[2] + 1.0
        image[18][1] = atom_with_image_dir[1] - 1.0
        image[18][2] = atom_with_image_dir[2] + 1.0

        # change 3 cord (+x, +y, +z) (-x, -y, -z) (-x, +y, +z) (+x, -y, +z)
        #               (+x, +y, -z) (+x, -y, -z) (-x, +y, -z) (-x, -y, +z)
        image[19][0] = atom_with_image_dir[0] + 1.0
        image[19][1] = atom_with_image_dir[1] + 1.0
        image[19][2] = atom_with_image_dir[2] + 1.0

        image[20][0] = atom_with_image_dir[0] - 1.0
        image[20][1] = atom_with_image_dir[1] - 1.0
        image[20][2] = atom_with_image_dir[2] - 1.0

        image[21][0] = atom_with_image_dir[0] - 1.0
        image[21][1] = atom_with_image_dir[1] + 1.0
        image[21][2] = atom_with_image_dir[2] + 1.0

        image[22][0] = atom_with_image_dir[0] + 1.0
        image[22][1] = atom_with_image_dir[1] - 1.0
        image[22][2] = atom_with_image_dir[2] + 1.0

        image[23][0] = atom_with_image_dir[0] + 1.0
        image[23][1] = atom_with_image_dir[1] + 1.0
        image[23][2] = atom_with_image_dir[2] - 1.0

        image[24][0] = atom_with_image_dir[0] + 1.0
        image[24][1] = atom_with_image_dir[1] - 1.0
        image[24][2] = atom_with_image_dir[2] - 1.0

        image[25][0] = atom_with_image_dir[0] - 1.0
        image[25][1] = atom_with_image_dir[1] + 1.0
        image[25][2] = atom_with_image_dir[2] - 1.0

        image[26][0] = atom_with_image_dir[0] - 1.0
        image[26][1] = atom_with_image_dir[1] - 1.0
        image[26][2] = atom_with_image_dir[2] + 1.0
        return image

    def close_image_dis(self, single_atom):
        dis = np.zeros(27)
        image = AtomImage(self.atom_with_image, self.cell).image()
        image_car = direct_cartesian_transform(image, self.cell, 'DtoC')
        for d in range(27):
            dis[d] = distance(image_car[d], single_atom)
        close_image = np.amin(dis)
        return close_image
    
    # single atom being the first atom in the rigid body 
    def close_image_position(self, single_atom):
        dis = np.zeros(27)
        image = AtomImage(self.atom_with_image, self.cell).image()
        image_car = direct_cartesian_transform(image, self.cell, 'DtoC')
        for d in range(27):
            dis[d] = distance(image_car[d], single_atom)
        close_image_index = int(np.where(dis == np.amin(dis))[0][0])
        close_image_position = image_car[close_image_index]
        return close_image_position


def dis_check(atom_positions, cell, dis_limit):
    length = len(atom_positions)
    count2 = 0
    for i in range(length - 1):
        count1 = 0
        for j in range(i + 1, length):
            dis = AtomImage(atom_positions[i], cell).close_image_dis(atom_positions[j])
            if dis >= dis_limit:
                count1 += 1
            # else:
            #     print(f'dis_check distance between atom{i+1} and atom{j+1} is {round(dis, 8)}')
        if count1 == length - i - 1:
            count2 += 1
    if count2 == length - 1:
        # print(f'count2:{count2}')
        # print('yes')
        return 'yes'
    else:
        # print('no')
        return 'no'
    # print(f'count2:{count2}')
    # print('-----------------------')


class WriteLog:
    def __init__(self, output_path):
        self.output_path = output_path

    def write_script_count(self):
        log_path = self.output_path + '/script_count'
        os.chdir(self.output_path)
        script1 = os.path.isfile(log_path)
        script = str(script1)
        if script == 'False':
            with open('script_count', 'w') as scrip:
                scrip.write('0')
        with open('script_count', 'r') as scrip_read:
            line_script = scrip_read.readlines()
            script_count = float(line_script[0].split()[0])
        script_count = script_count + 1
        with open('script_count', 'w') as scrip_write:
            scrip_write.truncate(0)
            scrip_write.write(f'{script_count}')
        return script_count
    
    def write_big_loop(self, script_count, internal_circulation):
        log_path = self.output_path + '/big_loop'
        os.chdir(self.output_path)
        script1 = os.path.isfile(log_path)
        script = str(script1)
        if script == 'False':
            with open('big_loop', 'w') as scrip:
                scrip.write('0')
        print(f'script_count:{script_count}')
        print(f'internal_circulation:{internal_circulation}')
        with open('big_loop', 'r') as scrip_read:
            line_script = scrip_read.readlines()
            big_loop = float(line_script[0].split()[0])
        if script_count % internal_circulation == 0:
            big_loop = big_loop + 1
            with open('big_loop', 'w') as scrip_write:
                scrip_write.truncate(0)
                scrip_write.write(f'{big_loop}')
        return big_loop
    
    def write_big_loop_plus(self):
        log_path = self.output_path + '/big_loop'
        os.chdir(self.output_path)
        script1 = os.path.isfile(log_path)
        script = str(script1)
        if script == 'False':
            with open('big_loop', 'w') as scrip:
                scrip.write('0')
        with open('big_loop', 'r') as scrip_read:
            line_script = scrip_read.readlines()
            big_loop = float(line_script[0].split()[0])
        big_loop = big_loop + 1
        with open('big_loop', 'w') as scrip_write:
            scrip_write.truncate(0)
            scrip_write.write(f'{big_loop}')
        return big_loop


def read_energy(keyword_file):
    keyword_text = 'sigma'
    with open(keyword_file, 'r') as outcar:
        lines = outcar.readlines()
        for line in reversed(lines):
            if keyword_text in line:
                energy = float(line.split()[6])
                break
    return energy


def metropolis_group_single(algorithm, energy, ref_energy, temp, metro_constant):
    accept = 0
    if algorithm == 'default':
        x = 1.0 / (1.0 + math.exp(((energy - ref_energy) * metro_constant)/ temp))
        random_number = random.uniform(0, 1)
        if x > 0.5:
            accept = 1
        elif x == 1:
            accept = 2
        elif x > random_number:
            accept = 1
        else:
            accept = 0
    if algorithm == 'test':
        x = math.exp((ref_energy - energy) / temp)
        random_number = random.uniform(0, 1)
        if x > 1:
            accept = 1
        elif x == 1:
            accept = 2
        elif x > random_number:
            accept = 1
        else:
            accept = 0
    return accept


def setup_file_copy(my_cwd, sym_no, case):
    os.chdir(my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/input_file')
    shutil.copy('INCAR', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/trail')
    shutil.copy('POTCAR', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/trail')
    shutil.copy('KPOINTS', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/trail')
    shutil.copy('vasp.sub', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/trail')
    shutil.copy('INCAR_relax', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/relax' + '/INCAR')
    shutil.copy('POTCAR', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/relax')
    shutil.copy('KPOINTS_relax', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/relax' + '/KPOINTS')
    shutil.copy('vasp.sub', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/relax')
    shutil.copy('INPUT', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('vasp2.sub', my_cwd + f'/SymTest/sym{sym_no}/{case}')

    os.chdir(my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/log_file')
    shutil.copy('POSCAR', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/trail')
    shutil.copy('POSCAR', my_cwd + f'/SymTest/sym{sym_no}/{case}' + '/relax')
    
    os.chdir(my_cwd)
    shutil.copy('functions.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('interprocess.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('metropolis.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('rigid.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('symdata.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('structureupdate.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('readinginput.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('randommove.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('mc.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('firstrelax.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')
    shutil.copy('forcejump.py', my_cwd + f'/SymTest/sym{sym_no}/{case}')