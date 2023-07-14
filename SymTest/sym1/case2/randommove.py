# randommove > rigid > function
from rigid import *

class AtomRandomMove:
    def __init__(self, initial_position, cell):
       self.initial_position = initial_position
       self.cell = cell
    
    def symmetry_restricted(self, rigid_or_single, algorithm, temp, wyckoff_position, *args):
        if rigid_or_single == 'rigid':
            rigid_type = rigid_select(args[0])
            center = rigid_type(1, 1.0).center(self.initial_position)
            atoms_move = np.array([0.0, 0.0, 0.0])
            center_shift = center.copy()
            if wyckoff_position[0] != 'x':
                atoms_move[0] = 0
            else:
                atoms_move[0] = step_algorithm(algorithm, temp)
            if wyckoff_position[1] != 'y':
                atoms_move[1] = 0
            else:
                atoms_move[1] = step_algorithm(algorithm, temp)
            if wyckoff_position[2] != 'z':
                atoms_move[2] = 0
            else:
                atoms_move[2] = step_algorithm(algorithm, temp)
            center_shift += atoms_move
            for i in range(3):
                center_shift[i] = round(center_shift[i], 8)
            final_position1 = self.initial_position + atoms_move
            final_position = move_atoms_into_box(final_position1, self.cell)
        
        elif rigid_or_single == 'single':
            atoms_move = np.array([.1, .1, .1])
            atoms_move_single = np.array([0.0, 0.0, 0.0])
            shift = self.initial_position.copy()
            if wyckoff_position[0] != 'x':
                atoms_move_single[0] = 0
            else:
                atoms_move_single[0] = step_algorithm(algorithm, temp)
            if wyckoff_position[1] != 'y':
                atoms_move_single[1] = 0
            else:
                atoms_move_single[1] = step_algorithm(algorithm, temp)
            if wyckoff_position[2] != 'z':
                atoms_move_single[2] = 0
            else:
                atoms_move_single[2] = step_algorithm(algorithm, temp)
            shift += atoms_move_single
            for j in range(3):
                shift[0][j] = round(shift[0][j], 8)
            atoms_move = atoms_move_single.copy()
            final_position1 = self.initial_position + atoms_move
            final_position = move_atoms_into_box(final_position1, self.cell)
               
        return final_position
    

# class cell_shape:
#     def __init__(self, initial_position, cell):
#        self.initial_position = initial_position
#        self.cell = cell