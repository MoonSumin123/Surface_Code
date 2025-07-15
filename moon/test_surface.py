import numpy as np
from simulation import *

class surface:
    def __init__(self):
        '''
        data_qubits = []
        X_stabilizers = []
        Z_stabilizers = []
        '''
        self.data_qubits = []
        self.X_stabilizers = []
        self.Z_stabilizers = []

    def add_data(self, q):
        self.data_qubits.append(q)

    def add_x(self, x):
        self.X_stabilizers.append(x)

    def add_z(self, z):
        self.Z_stabilizers.append(z)

class data_qubit:
    def __init__(self, x_boundary: bool = False, z_boundary: bool = False, t_x_boundary: bool = False, t_z_boundary: bool = False):
        """
        X_stabilizers = []
        Z_stabilizers = []
        is_x_boundary, is_z_boundary = whether this qubit lies on X/Z boundary
        """
        self.X_stabilizers = []
        self.Z_stabilizers = []
        self.is_x_boundary = x_boundary
        self.is_z_boundary = z_boundary
        self.value = 0

        self.test_x_b = t_x_boundary
        self.test_z_b = t_z_boundary

    def add_x(self, x):
        self.X_stabilizers.append(x)

    def add_z(self, z):
        self.Z_stabilizers.append(z)

class stabilizer:
    def __init__(self, boundary = False):
        '''
        data_qubits = []
        is_boundary = bool (True: bondary, False: None)
        index = integer of # of index of surface model's list
        '''
        self.data_qubits = []
        self.is_boundary = boundary
        self.index = 0
        self.value = 0

    def add_data(self, q):
        self.data_qubits.append(q)

    def set_ind(self, i):
        self.index = i

def build_rotated_surface_code(d):
    surf = surface()
    size = 2 * d - 1

    data_qubits = [[None for _ in range(size)] for _ in range(size)]

    # 1. 데이터 큐비트 생성 및 등록
    for i in range(size):
        for j in range(size):
            if (i + j) % 2 == 0:
                # 경계 조건 판별
                x_b = (i == 0 or i == size - 1)
                z_b = (j == 0 or j == size - 1)
                t_x_b = (i == 0)
                t_z_b = (j == 0)
                q = data_qubit(x_boundary=x_b, z_boundary=z_b, t_x_boundary=t_x_b, t_z_boundary=t_z_b)

                data_qubits[i][j] = q
                surf.add_data(q)

    # 2. Stabilizer 생성 및 연결
    for i in range(size):
        for j in range(size):
            neighbors = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
            qubit_neighbors = [
                (ni, nj)
                for ni, nj in neighbors
                if 0 <= ni < size and 0 <= nj < size and data_qubits[ni][nj] is not None
            ]

            # X stabilizer
            if i % 2 == 0 and j % 2 == 1:
                is_boundary = 1 if j == 1 or j == size-2 else 0
                stab = stabilizer(boundary=is_boundary)
                stab.set_ind(len(surf.X_stabilizers))
                for ni, nj in qubit_neighbors:
                    q = data_qubits[ni][nj]
                    stab.add_data(q)
                    q.add_x(stab)
                surf.add_x(stab)

            # Z stabilizer
            if i % 2 == 1 and j % 2 == 0:
                is_boundary = 1 if i == 1 or i == size-2 else 0
                stab = stabilizer(boundary=is_boundary)
                stab.set_ind(len(surf.Z_stabilizers))
                for ni, nj in qubit_neighbors:
                    q = data_qubits[ni][nj]
                    stab.add_data(q)
                    q.add_z(stab)
                surf.add_z(stab)

    return surf


list_of_distance = [3, 5, 7, 9, 11]
list_of_probability = [1e-2, 2e-2, 4e-2, 5e-2, 6e-2, 7e-2, 1e-1]
iteration = 3000

surf = build_rotated_surface_code(5)

result_total = []   # [[d1, [[p1, rate11], [p2,rate12], ...]],
                        #  [d2, [[p1, rate21], [p2,rate22], ...]], ... ]
for distance in list_of_distance:
    result_distance = []
    for probability in list_of_probability:
        logical_error = 0
        for i in range(iteration):
            target_surface = build_rotated_surface_code(distance)
            logical_error += one_cycle_simulation(target_surface, probability)
        rate_logical_error = logical_error / iteration
        result_distance.append([probability, rate_logical_error])
    result_total.append([distance, result_distance])

make_chart(result_total)