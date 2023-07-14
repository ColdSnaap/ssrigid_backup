from ase.spacegroup import crystal
from readinginput import *
from randommove import *

root_path = os.getcwd()
relax_path = root_path + '/relax'
trail_path = root_path + '/trail'
output_path = root_path + '/output'
log_file_path = root_path + '/log_file'
good_ones_relaxed_path = output_path + '/good_ones_relaxed'

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

# save good structures
isdir1 = os.path.isdir(good_ones_relaxed_path)
isdir = str(isdir1)
if isdir == 'False':
    os.chdir(output_path)
    os.mkdir('good_ones_relaxed')
os.chdir(relax_path)
energy1 = read_energy('OUTCAR')
energy = round(float(energy1), 4)
list_file = os.listdir(good_ones_relaxed_path)
number_files = len(list_file)
if number_files < 11:
    shutil.copytree(relax_path, good_ones_relaxed_path + '/' + str(-energy))
elif number_files == 11:
    file = []
    for i in range(number_files):
        file.append(float(list_file[i]))
    file_sorted = sorted(file)
    shutil.rmtree(good_ones_relaxed_path + f'/{str(file_sorted[0])}')

# relaxed vs script_count
os.chdir(output_path)
with open('script_count', 'r') as scrip_read:
    line_script = scrip_read.readlines()
    script_count = float(line_script[0].split()[0])
with open('energy_relax', 'a') as en_relax:
    en_relax.write(f'{script_count}  {energy}' + '\n')

os.chdir(root_path)