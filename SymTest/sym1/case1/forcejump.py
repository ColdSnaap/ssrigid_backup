# mc without log write

from ase.spacegroup import crystal
from readinginput import *
from randommove import *
from symdata import *

root_path = os.getcwd()
relax_path = root_path + '/relax'
trail_path = root_path + '/trail'
output_path = root_path + '/output'
log_file_path = root_path + '/log_file'
good_ones_relaxed_path = output_path + '/good_ones_relaxed'

"""
mc.py
"""
os.chdir(log_file_path)
get_no = RigidSingleNo()
rigid_no = get_no.rigid_no
single_no = get_no.single_no
WyckoffRigid = read_file('sym_information', 'WyckoffRigid')
WyckoffSingle = read_file('sym_information', 'WyckoffSingle')
os.chdir(root_path)

# log
os.chdir('output')
with open('big_loop', 'r') as scrip_read:
    line_script = scrip_read.readlines()
big_loop = float(line_script[0].split()[0])

# read POSCAR
os.chdir(trail_path)
poscar = ReadPOSCAR('POSCAR')
cell_poscar = poscar.lattice
poscar_position = poscar.positions

# current temp
temp = temp_algorithm(temp_update, temp_star, big_loop)

# generate a rigid class for center method
rigid_type_class = rigid_select(rigid_type)(1.0, 1.0)

close_dis_check = 'no'

print('forced jump ...')
while close_dis_check != 'yes':
    os.chdir(log_file_path)
    contcar_to_rigid(poscar_position)
    for n in range(rigid_no):
        wyckoff_position_rigid = sg[f'sg_{sym_no}'][WyckoffRigid[n]][2]
        wyckoff_element_rigid = sg[f'sg_{sym_no}'][WyckoffRigid[n]][1]
        rigid_array = rigid_to_array(f'rigid_{n+1}')
        rigid_after_move1 = AtomRandomMove(rigid_array, cell_poscar).symmetry_restricted('rigid', step_update, temp, wyckoff_position_rigid, rigid_type)
        # image select
        rigid_after_move = rigid_after_move1.copy()
        for i in range(1, len(rigid_after_move)):
            image_select = AtomImage(rigid_after_move1[i], cell_poscar).close_image_position(rigid_after_move1[0])
            rigid_after_move[i] = image_select
        
        # rotation
        rigid_center = rigid_type_class.center(rigid_after_move)
        rigid_after_rotation = rigid_rotation(rigid_after_move, rigid_center).random_rotation(sym_no, wyckoff_element_rigid, 10, cell_poscar)
        rigid_array_final = rigid_after_rotation
        array_to_rigid(rigid_array_final, f'rigid_{n+1}')
        write_rigid_trans(f'rigid_{n+1}', sym_no, cell_poscar)

    for m in range(single_no):
        wyckoff_position_single = sg[f'sg_{sym_no}'][WyckoffSingle[m]][2]
        single_array = rigid_to_array(f'single_{m+1}')
        single_after_move = AtomRandomMove(single_array, cell_poscar).symmetry_restricted('single', step_update, temp, wyckoff_position_single)
        array_to_rigid(single_after_move, f'single_{m+1}')
        write_rigid_trans(f'single_{m+1}', sym_no, cell_poscar)

    poscar_dir = write_poscar_dic()
    dir_to_poscar(poscar_dir, cell_poscar)
    # atom distance check
    poscar_read = ReadPOSCAR('POSCAR')
    atom_positions = poscar_read.positions   
    close_dis_check = dis_check(atom_positions, cell_poscar, dis_limit)

# copy poscar
shutil.copy('POSCAR', trail_path + '/POSCAR')
print('jump completed')
