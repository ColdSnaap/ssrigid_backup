# initial_cell function only has cube and hex, universal vol 
"""
Class SymCases has the following directory format
The example is with tetrahedron and symmetry number 15)
{'1': [5, [['s1', 's6'], ['s2', 's6'], ['s3', 's6'], ['s4', 's6'], ['s5', 's6']]],
'2': [5, [['s1', 's2', 's3'], ['s1', 's2', 's4'], ['s1', 's2', 's5'], ['s1', 's3', 's4'], ['s1', 's3', 's5'],
['s1', 's4', 's5'], ['s2', 's3', 's4'], ['s2', 's3', 's5'], ['s2', 's4', 's5'], ['s3', 's4', 's5']]]}

'1' is the case index means the 1st case is with the rigid body at 5th Wyckoff position
['s1', 's6'] means single atoms are at s1 and s6

Add new rigid body: Add new rigid body information in rigid.py and rigid_generate function in functions.py
"""

from ase.spacegroup import crystal
from readinginput import *
from rigid import *
from casecount import *

case = SymCases(rigid_type, ratio, 2)
case_dir = case.case_list(sym_no)

print(compound_single_atom)

rigid = rigid_select(rigid_type)
# trail_len = 0 means there is no possible combination of atoms
my_cwd = os.getcwd()
# make the test file where all the cases will be put in
path = my_cwd + '/SymTest'
isdir1 = os.path.isdir(path)
isdir = str(isdir1)
if isdir == 'False':
    os.mkdir('SymTest')

os.chdir(my_cwd + '/SymTest')
path_case = my_cwd + '/SymTest/' + f'sym{sym_no}'
isdir1 = os.path.isdir(path_case)
isdir = str(isdir1)
if isdir == 'False':
    os.mkdir(f'sym{sym_no}')
else:
    shutil.rmtree(f'{os.getcwd()}/sym{sym_no}', ignore_errors=True)
    # os.rmdir(f'sym{sym_no}')
    os.mkdir(f'sym{sym_no}')
my_cwd_cases = os.getcwd()

# generate single atoms
# total number of combination

os.chdir(my_cwd_cases + f'/sym{sym_no}')

for key in case_dir.keys():
    os.mkdir(key)
    os.chdir(os.getcwd() + f'/{key}')
    os.mkdir('trail')
    os.mkdir('relax')

    # write sym_information
    os.mkdir(f'log_file')
    os.chdir(os.getcwd() + '/log_file')
    sym_information(sym_no, case_dir, key)

    dis_transform_return = 'no'
    density_adjust = density

    # total atom number in rigid
    WyckoffRigid = read_file('sym_information', 'WyckoffRigid')
    WyckoffSingle = read_file('sym_information', 'WyckoffSingle')
    def total_atom_number(sym_no, wyckoff_rigid, wyckoff_single):
        atom_in_rigid = len(rigid(1, 1).atoms_position())
        number = 0
        for i in wyckoff_rigid:
            muti = sg[f'sg_{sym_no}'][i][0]
            number += atom_in_rigid * muti
        for j in wyckoff_single:
            muti = sg[f'sg_{sym_no}'][j][0]
            number += muti
        return number

    total_atom = total_atom_number(sym_no, WyckoffRigid, WyckoffSingle)

    dis_transform_return == 'no'

    while dis_transform_return != 'yes':
        density_adjust = density_adjust - 0.002
        count = 0
        while dis_transform_return != 'yes':
            # sim box
            sim_box = initial_cell(box_type, total_atom, density_adjust)
            # print(f'sim_box:{sim_box}')
            # generate rigid body
            
            n_rigid = 1
            for i in WyckoffRigid:
                wyckoff_positions = sg[f'sg_{sym_no}'][i][2]
                wyckoff_element = sg[f'sg_{sym_no}'][i][1]
                rigid_generate_inter = rigid(sym_no, bond).atoms_position() @ rigid(sym_no, bond).rotaiton_matrix(wyckoff_element) @ setup_random_rotation(sym_no, wyckoff_element, 360.0, sim_box)
                rigid_generate = move_rigid_to_wyckoff(rigid_generate_inter, wyckoff_positions, rigid_type, sim_box)
                with open(f'rigid_{n_rigid}', 'w') as rigid_write:
                    for j in compound_rigid:
                        rigid_write.write(f'{j} ')
                    rigid_write.write('\n')
                    for j in range(len(rigid_generate)):
                        rigid_write.write(f'{rigid_generate[j][0]} {rigid_generate[j][1]} {rigid_generate[j][2]}' + '\n')
                n_rigid += 1
            
            # generate single atoms
            n_single = 1
            for i in WyckoffSingle:
                wyckoff_positions_single = sg[f'sg_{sym_no}'][i][2]
                single_initial = generate_single_atoms(wyckoff_positions_single, sim_box)
                with open(f'single_{n_single}', 'w') as single_write:
                    for j in compound_single_atom:
                        single_write.write(f'{j} ')
                    single_write.write('\n')
                    single_write.write(f'{single_initial[0]} {single_initial[1]} {single_initial[2]}' + '\n')
                n_single += 1
        
            for i in range(1, n_rigid):
                write_rigid_trans(f'rigid_{i}', sym_no, sim_box)
            for i in range(1, n_single):
                write_rigid_trans(f'single_{i}', sym_no, sim_box)
            
            poscar_dir = write_poscar_dic()
            dir_to_poscar(poscar_dir, sim_box)

            # atom distance check
            poscar_read = ReadPOSCAR('POSCAR')
            atom_positions = poscar_read.positions
            
            dis_transform_return = dis_check(atom_positions, sim_box, dis_limit)
            count += 1
            if count == 300:
                break
        
        print(f'count:{count}')
        print(f'density_adjust:{density_adjust}')

    cell_final = ReadPOSCAR('POSCAR').lattice
    for i in range(1, n_rigid):
        write_atom_index(f'rigid_{i}', cell_final)
    for i in range(1, n_single):
        write_atom_index(f'single_{i}', cell_final)

    shutil.copytree(my_cwd + '/input_file', my_cwd_cases + f'/sym{sym_no}/{key}/input_file')


    # copy file to individual trail
    setup_file_copy(my_cwd, sym_no, key)

    os.chdir(my_cwd_cases + f'/sym{sym_no}')

# setup all input file 
# write poscar and atom position before transform
os.chdir(my_cwd)

