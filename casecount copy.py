from rigid import *
from symdata import *
from itertools import combinations, combinations_with_replacement

# PS4LiCa

def rigid_wyckoff_sym(rigid_type, sym_no):
        if rigid_type == 'Tetrahedron':
            wyckoff_sym = Tetrahedron(sym_no, 2.0).wyckoff_sym()
        if rigid_type == 'GeSe':
            wyckoff_sym = GeSe(sym_no, 2.0).wyckoff_sym()
        return wyckoff_sym


def total_atom_number_from_wyckoff(sym_no, wyckoff_list):
    count = 0
    for i in wyckoff_list:
        muti = sg[f'sg_{sym_no}'][i][0]
        count = count + muti
    return count


def uniq_site(sym_no, site):
    count = 0
    for i in sg[f'sg_{sym_no}'][site][2]:
        if i != 'x' and i != 'y' and i != 'z':
            count += 1
    if count == 3:
        return 'yes'
    else:
        return 'no'


rigid_max = 3
sym_no = 1
ratio = [1,1,1]
rigid_type = 'Tetrahedron'
# ratio = 1
# print(rigid_wyckoff_sym(rigid_type, 31))
for sym_no in range(1, 74):
    rigid_list = []
    single1_list = []
    single2_list = []

    sg_len = len(sg[f'sg_{sym_no}'])
    sg_list = []

    for i in range(1, sg_len+1):
        sg_list.append(f's{i}')
    sg_list_fi = []
    for i in sg_list:
        if sg[f'sg_{sym_no}'][i][1] in rigid_wyckoff_sym(rigid_type, sym_no):
            sg_list_fi.append(i)
    # print(sg_list_fi)
    sg_len_fi = len(sg_list_fi)

    for i in range(1, rigid_max+1):
        com_rigid= list(combinations_with_replacement(sg_list_fi, i))
        for j in com_rigid:
            rigid_list.append(j)
    # rigid_list.append(('s1','s1','s2','s2'))

    # get rid of single site and the ones rigid body number exceed the max number allowed
    # rigid_list_copy1 = no sighle site
    # rigid_list_copy2 = no single site and no exceed max rigid body number
    rigid_list_copy1 = rigid_list.copy()
    # print(len(com_rigid[0]))
    for y in rigid_list:
        dup = [x for i, x in enumerate(y) if i != y.index(x)]
        if len(dup) != 0:
            for i in dup:
                if uniq_site(sym_no, i) == 'yes':
                    if y in rigid_list_copy1:
                        rigid_list_copy1.remove(y)

    # print(rigid_list_copy1)
    rigid_list_copy2 = rigid_list_copy1.copy()
    for i in rigid_list_copy1:
        if total_atom_number_from_wyckoff(sym_no, i) > rigid_max:
            rigid_list_copy2.remove(i)

    rigid_fi = rigid_list_copy2.copy()
    # print(f'rigid_fi:{rigid_fi}')
    # # -----------------------------------------------------------------
    count3 = 0
    # single2
    # print('-------------------------------------')
    for rigid_i in rigid_fi:
        count4 = 0
        if len(rigid_i) == 1:
            muti = total_atom_number_from_wyckoff(sym_no, [rigid_i[0]])
        else:
            muti = total_atom_number_from_wyckoff(sym_no, rigid_i)
        single1_max = muti*ratio[1]

        for k in range(1, single1_max+1):
            com_single= list(combinations_with_replacement(sg_list, k))
            for j in com_single:
                single1_list.append(j)
        # print(rigid_list)
        single1_list_copy1 = single1_list.copy()
        # print(len(com_rigid[0]))
        for y in single1_list:
            dup = [x for i, x in enumerate(y) if i != y.index(x)]
            if len(dup) != 0:
                for i in dup:
                    if uniq_site(sym_no, i) == 'yes':
                        if y in single1_list_copy1:
                            single1_list_copy1.remove(y)

        single1_list_copy2 = single1_list_copy1.copy()
        # print(rigid_fi)
        # print(single_list_copy2)
        # print(single_max)
        for h in single1_list_copy1:
            if total_atom_number_from_wyckoff(sym_no, h) != single1_max:
                single1_list_copy2.remove(h)

        single1_fi = single1_list_copy2.copy()
        # a Z      print(f'rigid:{rigid_i} single1:{single1_fi}')
        if len(ratio) == 3:
            single2_max = muti*ratio[2]

            for k in range(1, single2_max+1):
                com_single= list(combinations_with_replacement(sg_list, k))
                for j in com_single:
                    single2_list.append(j)
            # print(rigid_list)
            single2_list_copy1 = single2_list.copy()
            # print(len(com_rigid[0]))
            for y in single2_list:
                dup = [x for i, x in enumerate(y) if i != y.index(x)]
                if len(dup) != 0:
                    for i in dup:
                        if uniq_site(sym_no, i) == 'yes':
                            if y in single2_list_copy1:
                                single2_list_copy1.remove(y)

            single2_list_copy2 = single2_list_copy1.copy()
            for h in single2_list_copy1:
                if total_atom_number_from_wyckoff(sym_no, h) != single2_max:
                    single2_list_copy2.remove(h)

            single2_fi = single2_list_copy2.copy()
            count4 += len(single1_fi)*len(single2_fi)
            count3 += count4
            # print(f'single2:{single2_fi}')
            # print('----------------------------')

    # different type of single atom at uniq position
        


    #     count3 += len(single_fi)*len(single_fi)

    print(f'{sym_no}  {count3}')








# same initial rigibody and single atom can only be at one wyckoff position 
# class SymCases:
#     def __init__(self, rigid_type, ratio, rigid_number_limit):
#         self.rigid_type = rigid_type
#         self.ratio = ratio
#         self.rigid_no_limit = rigid_number_limit
    
#     # rigid body can't be at the same wyckoff position except for general position
#     # works on single type of rigid bodies, single type of single atom
#     def case_list(self, sym_no):
#         case_dic = {}
#         rigid_list = []
#         sym_len = len(sg[f'sg_{sym_no}'])
#         # the max number allowed at general position number_limit / general_muti
#         general_muti = sg[f'sg_{sym_no}'][f's{sym_len}'][0]
#         max_number = int(self.rigid_no_limit / general_muti)
#         # conbination of 1
#         for i in range(sym_len - 1):
#             i_sym = sg[f'sg_{sym_no}'][f's{i + 1}'][1]
#             if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
#                 rigid_list.append([f's{i + 1}'])
#         # conbination of 2 to sym_len - 1
#         # print(rigid_list)
#         special_wyckoff_list = []
        
#         # symmetry
#         for i in range(sym_len - 1):
#             i_sym = sg[f'sg_{sym_no}'][f's{i + 1}'][1]
#             if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
#                 special_wyckoff_list.append(f's{i + 1}')
#         # print(f'special_wyckoff_list:{special_wyckoff_list}')
#         # for i in range(sym_len - 1):
#         #     special_wyckoff_list.append(f's{i + 1}')
        
#         sym_len_special = len(special_wyckoff_list)
#         # print(f'special_wyckoff_list:{special_wyckoff_list}')
#         # comb of sym_len_special
#         # print(f'rigid_list:{rigid_list}')
#         for i in range(2, sym_len_special+1):
#             special_combine = list(combinations(special_wyckoff_list, i))
#             for j in special_combine:
#                 rigid_list.append(list(j))
#         # print(f'rigid_list:{rigid_list}')
#         # print(special_combine)
#         # get rid of the ones exceed the max_muti
#         # print(rigid_list)
#         # print(len(rigid_list))
#         rigid_list_copy = rigid_list.copy()

#         for i in rigid_list_copy:
#             check1 = total_atom_number_from_wyckoff(sym_no, i)
#             if check1 > self.rigid_no_limit:
#                 rigid_list.remove(i)
#         # take general position into consideration
#         if max_number == 0:
#             list_final = rigid_list.copy()
#         else:
#             rigid_list.append([f's{sym_len}'])
#             rigid_list_inter = rigid_list.copy()
#             # print(f'rigid_list:{rigid_list}')
#             # print(f'max_no:{max_number}')
#             for i in range(max_number):
#                 for j in rigid_list_inter:
#                     list_inter = j.copy()
#                     # print(f'i:{i}')
#                     # print(f'list_inter:{list_inter}')
#                     for k in range(i+1):
#                         list_inter.append(f's{sym_len}')
#                     # print(f'append:{list_inter}')
#                     check2 = total_atom_number_from_wyckoff(sym_no, list_inter)
#                     # print(f'check2:{check2}')
#                     # print(f'rigid_no_limit:{self.rigid_no_limit}')
#                     if check2 <= self.rigid_no_limit:
#                         rigid_list.append(list_inter)
#                     # print(f'rigid_list:{rigid_list}')
#                     # print(f'-------------------------------------')
#         list_final = rigid_list.copy()
#         # print(f'list_final:{list_final}')          
#         # check symmetry
#         # sym_check = list_final.copy()
#         # for i in sym_check:
#         #     i_len = len(i)
#         #     for j in range(i_len):
#         #         i_sym = sg[f'sg_{sym_no}'][i[j]][1]
#         #         if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
#         #             check_sym = 1
#         #         else:
#         #             check_sym = 0
#         #             break
#         #     if check_sym == 1:
#         #         pass
#         #     elif check_sym == 0:
#         #         list_final.remove(i)
#         # print(f'sym_check:{list_final}')

#         # get single list
#         if len(list_final) > 100:
#             case_dic['case1'] = 'TBD'
#         else:
#             list_final_inter = list_final.copy()
#             for i in list_final_inter:
#                 # if i != ['s1']:
#                 #     break
#                 single_muti = self.ratio * total_atom_number_from_wyckoff(sym_no, i)
#                 max_number_single = int(single_muti / general_muti)
#                 # conbination of 1
#                 single_list = []
#                 for j in range(sym_len - 1):
#                     single_list.append([f's{j + 1}'])
#                 # conbination of 2 to min(sym_len - 1, single_muti)
#                 for j in range(2, min(sym_len - 1, single_muti) + 1):
#                     special_combine_single = list(combinations(special_wyckoff_list, j))
#                     for k in special_combine_single:
#                         single_list.append(list(k))
#                 # print(f'special_wyckoff_list:{special_wyckoff_list}')
#                 # print(f'rigid:{i}')
#                 # print(f'single_muti:{single_muti}')
#                 # print(f'single_list:{single_list}')
#                 # print('-------------------------------')
#                 # get rid of the ones exceed the single_muti
#                 single_list_copy = single_list.copy()
#                 for j in single_list_copy:
#                     check3 = total_atom_number_from_wyckoff(sym_no, j)
#                     if check3 > single_muti:
#                         single_list.remove(j)
#                     for k in j:
#                         for l in i:
#                             if sg[f'sg_{sym_no}'][k][2][0] != 'x' and \
#                                 sg[f'sg_{sym_no}'][k][2][1] != 'y' and \
#                                 sg[f'sg_{sym_no}'][k][2][2] != 'z' and k == l:
#                                 single_list.remove(j)
#                                 break
#                 # print(f'single_list:{single_list}')
#                 # take general position into consideration
#                 if max_number_single == 0:
#                     pass
#                 else:
#                     single_list.append([f's{sym_len}'])
#                     # print(f'single_list:{single_list}')
#                     # print('-------------------------')
#                     single_list_inter = single_list.copy()
#                     for j in range(max_number_single):
#                         for k in single_list_inter:
#                             list_inter2 = k.copy()
#                             # print(f'k:{k}')
#                             for n in range(j+1):
#                                 list_inter2.append(f's{sym_len}')
#                             check4 = total_atom_number_from_wyckoff(sym_no, list_inter2)
#                             # print(f'max_number_single:{max_number_single}')
#                             # print(f'j:{j}')
#                             # print(f'single_list_inter:{list_inter2}')
#                             # print(f'single_list:{single_list}')
#                             # print(f'single_muti:{single_muti}')
#                             # print('---------------------------------')
#                             if check4 <= single_muti:
#                                 single_list.append(list_inter2)
#                             # print(f'j:{j}')
#                             # print(f'single_list:{single_list}')
#                             # print('---------------------------------')
#                 # ratio
#                 list_inter3 = single_list.copy()
#                 for j in list_inter3:
#                     check5 = total_atom_number_from_wyckoff(sym_no, j)
#                     if check5 != single_muti:
#                         single_list.remove(j)
                
#                 list_final_single = single_list.copy()
                
#                 index_case_dic = len(case_dic)
#                 index_list = len(list_final_single)
#                 for j in range(index_list):
#                     case_dic[f'case{index_case_dic + j + 1}'] = {}
#                     case_dic[f'case{index_case_dic + j + 1}']['rigid'] = i
#                     case_dic[f'case{index_case_dic + j + 1}']['single'] = list_final_single[j]
#         return case_dic    
    
# def case_count(rigid_type, ratio, rigid_number_limit, low, high):
#     sym = SymCases(rigid_type, ratio, rigid_number_limit)
#     for i in range(low, high + 1):
#         sym_dic = sym.case_list(i)
#         if len(sym_dic) == 1 and sym_dic['case1'] == 'TBD':
#             count = 'TBD'
#         elif len(sym_dic) == 0:
#             count = 0
#         else:
#             count = len(sym_dic)
#         print(f'sym:{i}  case:{count}')