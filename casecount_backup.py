from rigid import *
from symdata import *
from itertools import combinations, combinations_with_replacement

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


# same initial rigibody and single atom can only be at one wyckoff position 
class SymCases:
    def __init__(self, rigid_type, ratio, rigid_number_limit):
        self.rigid_type = rigid_type
        self.ratio = ratio
        self.rigid_no_limit = rigid_number_limit
    
    # rigid body can't be at the same wyckoff position except for general position
    # works on single type of rigid bodies, single type of single atom
    def case_list(self, sym_no):
        case_dic = {}
        rigid_list = []
        sym_len = len(sg[f'sg_{sym_no}'])
        # the max number allowed at general position number_limit / general_muti
        general_muti = sg[f'sg_{sym_no}'][f's{sym_len}'][0]
        max_number = int(self.rigid_no_limit / general_muti)
        # conbination of 1
        for i in range(sym_len - 1):
            i_sym = sg[f'sg_{sym_no}'][f's{i + 1}'][1]
            if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
                rigid_list.append([f's{i + 1}'])
        # conbination of 2 to sym_len - 1
        # print(rigid_list)
        special_wyckoff_list = []
        
        # symmetry
        for i in range(sym_len - 1):
            i_sym = sg[f'sg_{sym_no}'][f's{i + 1}'][1]
            if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
                special_wyckoff_list.append(f's{i + 1}')
        # print(f'special_wyckoff_list:{special_wyckoff_list}')
        # for i in range(sym_len - 1):
        #     special_wyckoff_list.append(f's{i + 1}')
        
        sym_len_special = len(special_wyckoff_list)
        # print(f'special_wyckoff_list:{special_wyckoff_list}')
        # comb of sym_len_special
        # print(f'rigid_list:{rigid_list}')
        for i in range(2, sym_len_special+1):
            special_combine = list(combinations(special_wyckoff_list, i))
            for j in special_combine:
                rigid_list.append(list(j))
        # print(f'rigid_list:{rigid_list}')
        # print(special_combine)
        # get rid of the ones exceed the max_muti
        # print(rigid_list)
        # print(len(rigid_list))
        rigid_list_copy = rigid_list.copy()

        for i in rigid_list_copy:
            check1 = total_atom_number_from_wyckoff(sym_no, i)
            if check1 > self.rigid_no_limit:
                rigid_list.remove(i)
        # take general position into consideration
        if max_number == 0:
            list_final = rigid_list.copy()
        else:
            rigid_list.append([f's{sym_len}'])
            rigid_list_inter = rigid_list.copy()
            # print(f'rigid_list:{rigid_list}')
            # print(f'max_no:{max_number}')
            for i in range(max_number):
                for j in rigid_list_inter:
                    list_inter = j.copy()
                    # print(f'i:{i}')
                    # print(f'list_inter:{list_inter}')
                    for k in range(i+1):
                        list_inter.append(f's{sym_len}')
                    # print(f'append:{list_inter}')
                    check2 = total_atom_number_from_wyckoff(sym_no, list_inter)
                    # print(f'check2:{check2}')
                    # print(f'rigid_no_limit:{self.rigid_no_limit}')
                    if check2 <= self.rigid_no_limit:
                        rigid_list.append(list_inter)
                    # print(f'rigid_list:{rigid_list}')
                    # print(f'-------------------------------------')
        list_final = rigid_list.copy()
        # print(f'list_final:{list_final}')          
        # check symmetry
        # sym_check = list_final.copy()
        # for i in sym_check:
        #     i_len = len(i)
        #     for j in range(i_len):
        #         i_sym = sg[f'sg_{sym_no}'][i[j]][1]
        #         if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
        #             check_sym = 1
        #         else:
        #             check_sym = 0
        #             break
        #     if check_sym == 1:
        #         pass
        #     elif check_sym == 0:
        #         list_final.remove(i)
        # print(f'sym_check:{list_final}')

        # get single list
        if len(list_final) > 100:
            case_dic['case1'] = 'TBD'
        else:
            list_final_inter = list_final.copy()
            for i in list_final_inter:
                # if i != ['s1']:
                #     break
                single_muti = self.ratio * total_atom_number_from_wyckoff(sym_no, i)
                max_number_single = int(single_muti / general_muti)
                # conbination of 1
                single_list = []
                for j in range(sym_len - 1):
                    single_list.append([f's{j + 1}'])
                # conbination of 2 to min(sym_len - 1, single_muti)
                for j in range(2, min(sym_len - 1, single_muti) + 1):
                    special_combine_single = list(combinations(special_wyckoff_list, j))
                    for k in special_combine_single:
                        single_list.append(list(k))
                # print(f'special_wyckoff_list:{special_wyckoff_list}')
                # print(f'rigid:{i}')
                # print(f'single_muti:{single_muti}')
                # print(f'single_list:{single_list}')
                # print('-------------------------------')
                # get rid of the ones exceed the single_muti
                single_list_copy = single_list.copy()
                for j in single_list_copy:
                    check3 = total_atom_number_from_wyckoff(sym_no, j)
                    if check3 > single_muti:
                        single_list.remove(j)
                    for k in j:
                        for l in i:
                            if sg[f'sg_{sym_no}'][k][2][0] != 'x' and \
                                sg[f'sg_{sym_no}'][k][2][1] != 'y' and \
                                sg[f'sg_{sym_no}'][k][2][2] != 'z' and k == l and \
                                j in single_list:
                                single_list.remove(j)
                                break
                # print(f'single_list:{single_list}')
                # take general position into consideration
                if max_number_single == 0:
                    pass
                else:
                    single_list.append([f's{sym_len}'])
                    # print(f'single_list:{single_list}')
                    # print('-------------------------')
                    single_list_inter = single_list.copy()
                    for j in range(max_number_single):
                        for k in single_list_inter:
                            list_inter2 = k.copy()
                            # print(f'k:{k}')
                            for n in range(j+1):
                                list_inter2.append(f's{sym_len}')
                            check4 = total_atom_number_from_wyckoff(sym_no, list_inter2)
                            # print(f'max_number_single:{max_number_single}')
                            # print(f'j:{j}')
                            # print(f'single_list_inter:{list_inter2}')
                            # print(f'single_list:{single_list}')
                            # print(f'single_muti:{single_muti}')
                            # print('---------------------------------')
                            if check4 <= single_muti:
                                single_list.append(list_inter2)
                            # print(f'j:{j}')
                            # print(f'single_list:{single_list}')
                            # print('---------------------------------')
                # ratio
                list_inter3 = single_list.copy()
                for j in list_inter3:
                    check5 = total_atom_number_from_wyckoff(sym_no, j)
                    if check5 != single_muti:
                        single_list.remove(j)
                
                list_final_single = single_list.copy()
                
                index_case_dic = len(case_dic)
                index_list = len(list_final_single)
                for j in range(index_list):
                    case_dic[f'case{index_case_dic + j + 1}'] = {}
                    case_dic[f'case{index_case_dic + j + 1}']['rigid'] = i
                    case_dic[f'case{index_case_dic + j + 1}']['single'] = list_final_single[j]
        return case_dic    
    
def case_count(rigid_type, ratio, rigid_number_limit, low, high):
    sym = SymCases(rigid_type, ratio, rigid_number_limit)
    for i in range(low, high + 1):
        sym_dic = sym.case_list(i)
        if len(sym_dic) == 1 and sym_dic['case1'] == 'TBD':
            count = 'TBD'
        elif len(sym_dic) == 0:
            count = 0
        else:
            count = len(sym_dic)
        print(f'sym:{i}  case:{count}')