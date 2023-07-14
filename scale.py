# from casecount2 import *
import matplotlib.pyplot as plt
from readinginput import *
from itertools import combinations, combinations_with_replacement
from ase.spacegroup import crystal
from casecount import *
import numpy as np
from functions import *
from cellcheck import lattice_angle

sim_box = np.array([[10.695, -2.119, 0.594],
                [-1.692, 5.842, 1.004],
                [0.029, 1.677, 5.2001]])

# compound = ['P', 'S', 'S', 'S', 'S']

# sym_no = 31

# system_initial_tuple = [(0.0000000000000000, 0.3193858310057314, 0.3004191256249643),
#                         (0.2201941407319539, 0.1687908743112779, 0.1891720597598417),
#                         (0.7798058592740501, 0.1687908743112779, 0.1891720597598417),
#                         (0.0000000000000000, 0.6213979744118994, 0.1856207017628689),
#                         (0.0000000000000000, 0.3119861906253145, 0.6394921055933835)]

# system_transform = crystal(compound, system_initial_tuple, spacegroup=sym_no, cell=sim_box)

# print(system_transform.positions)

x = SymCases(rigid_type, ratio, 3)
# print(f'sym:{sym_no}')
print(f'rigid_type:{rigid_type}')
b = x.case_list(6)
for key in b.keys():
    print(f'{key}: {b[key]}')

# sym_information(sym_no, b, 3)
# root_path = os.getcwd()
# output_path = root_path + '/output'
# trail_path = root_path + '/SymTest/sym1/case2/trail'
# relax_path = root_path + '/SymTest/sym1/case2/relax'
# log_path = root_path + '/SymTest/sym31/case1/log_file'


# cell = [(6.9336127435, 0., 0.),
#         (0., 6.9336127435, 0.),
#         (0., 0., 6.9336127435)]

# atom_cor1 = np.array([6.684002684734, 3.4321383080325, 1.3312536467520002])
# atom_cor2 = direct_cartesian_transform(atom_cor1, cell, 'CtoD')
# atom_cor3 = [(atom_cor2[0], atom_cor2[1], atom_cor2[2])]
# # print(atom_cor2)
# # print(atom_cor3)
# # exit()

# x = crystal('Li', atom_cor3, spacegroup=6, cell=cell)
# p = x.positions

# print(p)


plt.plot([1, 2, 3, 4], [1, 4, 9, 10], 'ro')
plt.axis([0, 6, 0, 10])
plt.show()