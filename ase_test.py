# generate a rigid body at given wyckoff position 
from rigid import *
space_number = 118
wyckoff_position = [0, 0, 0]
rigid_type = rigid_select('Tetrahedron')
# rigid_composition = ['P', 'S', 'S', 'S', 'S']

cell = [(6.9336127435, 0., 0.),
        (0., 6.9336127435, 0.),
        (0., 0., 6.9336127435)]

rigid = rigid_type(space_number, 2.0)
rigid_cor = rigid.atoms_position()
rigid_to_wyckoff = move_rigid_to_wyckoff(rigid_cor, wyckoff_position, 'Tetrahedron', cell)
rotation = rigid.rotaiton_matrix(['-4', 'y', 'z'])
rigid_af_rotate = rigid_to_wyckoff @ rotation

print(rigid_af_rotate)
print('---------------')
for cor in rigid_af_rotate:
    cor2 = direct_cartesian_transform(cor, cell, 'CtoD')
    cor3 = [(cor2[0], cor2[1], cor2[2])]
    x = crystal('Li', cor3, spacegroup=space_number, cell=cell)
    p = x.positions
    print(p)