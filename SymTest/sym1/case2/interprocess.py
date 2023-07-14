from functions import *

root_path = os.getcwd()
trail_path = root_path +'/trail'
output_path = root_path +'/output'
good_ones_path = output_path + '/good_ones'

os.chdir(output_path)
isdir1 = os.path.isdir(good_ones_path)
isdir = str(isdir1)
if isdir == 'False':
    os.mkdir('good_ones')

os.chdir(trail_path)
energy1 = read_energy('OUTCAR')
energy = round(float(energy1), 4)

os.chdir(root_path)
list_file = os.listdir(good_ones_path)
number_files = len(list_file)
if number_files < 11:
    shutil.copytree('trail', good_ones_path + '/' + str(-energy))
elif number_files == 11:
    file = []
    for i in range(number_files):
        file.append(float(list_file[i]))
    file_sorted = sorted(file)
    shutil.rmtree(output_path + f'/good_ones/{str(file_sorted[0])}')

os.chdir(root_path)