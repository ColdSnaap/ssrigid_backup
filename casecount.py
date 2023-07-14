from rigid import *
from symdata import *
from itertools import combinations, combinations_with_replacement
import collections

def rigid_wyckoff_sym(rigid_type, sym_no):
        if rigid_type == 'Tetrahedron':
            wyckoff_sym = Tetrahedron(sym_no, 2.0).wyckoff_sym()
        if rigid_type == 'GeSe':
            wyckoff_sym = GeSe(sym_no, 2.0).wyckoff_sym()
        return wyckoff_sym


def total_atom_number_from_wyckoff(sym_no, wyckoff_list):
    count = 0
    # print(f'wyckoff_list:{wyckoff_list}')
    for i in wyckoff_list:
        # print(f'i:{i}')
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


class SymCases:
    def __init__(self, rigid_type, ratio, rigid_number_limit):
        self.rigid_type = rigid_type
        self.ratio = ratio
        self.rigid_no_limit = rigid_number_limit
    # ratio = [1, 3] not support [1, 1, 1] right now
    def case_list(self, sym_no):
        case_dic = {}
        sym_site_list = []
        rigid_site_list = []
        rigid_com_list = []
        rigid_com_list_final = []
        single_site_list = []
        single_site_list_final = []
        sym_len = len(sg[f'sg_{sym_no}'])

        for i in range(sym_len):
            sym_site_list.append([f's{i + 1}'])

        # print(f'sym_site_list:{sym_site_list}')

        for i in range(sym_len):
            i_sym = sg[f'sg_{sym_no}'][f's{i + 1}'][1]
            if i_sym in rigid_wyckoff_sym(rigid_type, sym_no):
                rigid_site_list.append([f's{i + 1}'])

        # print(f'sym_no:{sym_no}')
        # print(f'rigid_no_limit:{self.rigid_no_limit}')
        # print(f'rigid_site_list:{rigid_site_list}')

        for i in range(1, self.rigid_no_limit+1):
            com_rigid = list(combinations_with_replacement(rigid_site_list, i))
            for j in com_rigid:
                if i == 1:
                    rigid_com_list.append(j[0])
                else:
                     site_list = []
                     for k in j:
                        site_list.append(k[0])
                     rigid_com_list.append(site_list)
        # print(f'rigid_com_list:{rigid_com_list}')
        # find duplicated sites
        for i in rigid_com_list:
            uniq_site_check = 0
            dup = [item for item, count in collections.Counter(i).items() if count > 1]
            if len(dup) != 0:
                for n in dup:
                    if uniq_site(sym_no, n) == 'yes':
                        uniq_site_check = 1
                        break
            if total_atom_number_from_wyckoff(sym_no, i) <= self.rigid_no_limit and uniq_site_check == 0:
            # if total_atom_number_from_wyckoff(sym_no, i) <= self.rigid_no_limit:
                rigid_com_list_final.append(i)
        # print(f'rigid_com_list_final:{rigid_com_list_final}')
            # get single comb in i
        for i in rigid_com_list_final:
            single_ratio = self.ratio[1]
            single_amount = total_atom_number_from_wyckoff(sym_no, i) * single_ratio
            single_site_list1 = []
            single_site_list_final = []
            for m in range(1, single_amount+1):
                single_site_list = list(combinations_with_replacement(sym_site_list, m))
                for j in single_site_list:
                    if len(j) == 1:
                        single_site_list1.append(j[0])
                    else:
                        site_list = []
                        for k in j:
                            site_list.append(k[0])
                        single_site_list1.append(site_list)
            # print(f'single_site_list1:{single_site_list1}')
            for j in single_site_list1:
                # print(f'f:{j}')
                uniq_site_check = 0

                dup = [item for item, count in collections.Counter(j).items() if count > 1]
                # print(f'dup:{dup}')
                if len(dup) != 0:
                    for n in dup:
                        if uniq_site(sym_no, n) == 'yes':
                            uniq_site_check = 1
                            break
              
                if total_atom_number_from_wyckoff(sym_no, j) == single_amount and uniq_site_check == 0:
                    single_site_list_final.append(j)
            
            # print(f'i:{i}')
            # print(f'single_site_list_final:{single_site_list_final}')

            index_case_dic = len(case_dic)
            index_list = len(single_site_list_final)
            for j in range(index_list):
                case_dic[f'case{index_case_dic + j + 1}'] = {}
                case_dic[f'case{index_case_dic + j + 1}']['rigid'] = i
                case_dic[f'case{index_case_dic + j + 1}']['single'] = single_site_list_final[j]
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