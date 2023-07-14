from functions import *
from readinginput import *

root_path = os.getcwd()
trail_path = root_path +'/trail'
output_path = root_path +'/output'
log_path = output_path +'/log'
os.chdir(trail_path)

# get energy from outcar and outcar_ref
ref_energy = read_energy('OUTCAR_ref')
energy = read_energy('OUTCAR')
print('reference energy : ' + str(ref_energy))
print('energy : ' + str(energy))

os.chdir(output_path)
with open('big_loop', 'r') as scrip_read:
    line_script = scrip_read.readlines()
    big_loop = float(line_script[0].split()[0])

with open('metro_constant', 'r') as read:
    sd = read.readlines()
metro_constant = float(sd[0].split()[0])

# metroplis
current_temp = temp_algorithm(temp_update, temp_star, big_loop)
accept = metropolis_group_single('default', energy, ref_energy, current_temp, metro_constant)
print('current temp : {}'.format(current_temp))

os.chdir(trail_path)
if accept == 1:
    shutil.copy('POSCAR', 'POSCAR_ref')
    shutil.copy('OUTCAR', 'OUTCAR_ref')
    print('accept : yes')
else:
    print('accept : no')

# calculate accept rate
os.chdir(output_path)
isfile1 = os.path.isfile(log_path)
isfile = str(isfile1)
if isfile == 'False':
    with open('log', 'a') as log_count:
        for i in range(101):
            log_count.write('0' + '\n')
        log_count.write('1' + '\n')

with open('log', 'r') as log_read:
    lines = log_read.readlines()

loop_current = int(lines[100])
current_step = int(lines[101])

with open('log', 'w') as log_edit:
    for i, line in enumerate(lines):
        if accept == 1 and i == loop_current:
            log_edit.writelines('1\n')
        elif accept == 0 and i == loop_current:
            log_edit.writelines('0\n')
        else:
            log_edit.writelines(line)

with open('log', 'r') as log_read_2:
    lines_2 = log_read_2.readlines()

zero_count = 0
one_count = 0
print(f'current step {current_step}')
if current_step < 100:
    for i in range(current_step):
        if str(lines_2[i]) == '0\n':
            zero_count += 1
        elif str(lines_2[i]) == '1\n':
            one_count += 1
else:
    for i in range(100):
        if str(lines_2[i]) == '0\n':
            zero_count += 1
        elif str(lines_2[i]) == '1\n':
            one_count += 1

accept_rate = round(one_count / (zero_count + one_count), 5)

print(f'current acceptance rate : {accept_rate * 100}%')

if loop_current < 99:
    loop_current += 1
else:
    loop_current = 0

current_step += 1

# write log
with open('log', 'w') as log_edit_2:
    for i, line in enumerate(lines_2):
        if i == 100:
            log_edit_2.writelines(str(loop_current) + '\n')
        elif i == 101:
            log_edit_2.writelines(str(current_step) + '\n')
        else:
            log_edit_2.writelines(line)

# write energy
# if accept == 1:
with open('energy', 'a') as en:
    en.write(f'{current_step}  {energy}' + '\n')

# energy accept
if accept == 1:
    with open('energy_acc', 'a') as en_acc:
        en_acc.write(f'{current_step}  {energy}' + '\n')
elif accept == 0:
    with open('energy_rej', 'a') as en_rej:
        en_rej.write(f'{current_step}  {energy}' + '\n')

# debug
if debug == 1:
    debug_path = output_path + '/debug'
    isdir1 = os.path.isdir(debug_path)
    isdir = str(isdir1)
    if isdir == 'False':
        os.mkdir('debug')
    os.chdir(trail_path)
    shutil.copy('POSCAR', output_path + f'/debug/POSCAR_{int(current_step)}')

# # if the acceptance rate < 10% bigloop =+ 1
# big_loop_log = WriteLog(output_path)
# if accept_rate < 0.1:
#     big_loop_plus = big_loop_log.write_big_loop_plus()

os.chdir(root_path)