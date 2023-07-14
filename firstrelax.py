# structure update + mc
# get the energy from contcar
# calculation metro constant

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
structureupdate.py
"""
os.chdir(relax_path)
with open('CONTCAR', 'r') as contcar_read:
    contcar_line = contcar_read.readlines()
contcar = ReadPOSCAR('CONTCAR')
after_relax = contcar.positions
cell_after = contcar.lattice
after_relax_car = direct_cartesian_transform(after_relax, cell_after, 'DtoC')

# update
os.chdir(log_file_path)
contcar_to_rigid(after_relax_car)

# write CONTCAR from relax to POSCAR in trail
# change the dir coor from contcar to car in POSCAR
os.chdir(trail_path)
with open('POSCAR', 'w') as poscar_write:
    for i in range(7):
        poscar_write.write(contcar_line[i])
    poscar_write.write('Car\n')
    for i in range(len(after_relax_car)):
        x = round(after_relax_car[i][0], 8)
        y = round(after_relax_car[i][1], 8)
        z = round(after_relax_car[i][2], 8)
        poscar_write.write(str(x) + '  ' + str(y) + '  ' + str(z) + '\n')

shutil.copy('POSCAR', log_file_path + '/POSCAR')

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

# output file
output_path = root_path +'/output'
isfile1 = os.path.isdir(output_path)
isfile = str(isfile1)
if isfile == 'False':
    os.mkdir('output')

# read POSCAR
os.chdir(trail_path)
poscar = ReadPOSCAR('POSCAR')
cell_poscar = poscar.lattice
poscar_position = poscar.positions

# current temp
temp = temp_algorithm(temp_update, temp_star, 0.0)

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

# metro constant
os.chdir(relax_path)
e_relax = read_energy('OUTCAR')
metro_constant = round(-76.5/float(e_relax), 6)
# write metro constant
os.chdir(output_path)
with open('metro_constant', 'w') as cons:
    cons.write(str(metro_constant))