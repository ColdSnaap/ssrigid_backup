import numpy as np
atom_cor1 = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
atom_cor2 = np.array([[1, 2, 1],[1, 1, 1],[1, 1, 1]])

list = []
for j,k in enumerate(atom_cor1):
    list_inter = []
    for i in atom_cor1[j]:
        list_inter.append(i)
    list.append(list_inter)

for j, k in enumerate(atom_cor1):    
    print(f'{k[0]} {k[1]} {k[2]}')
