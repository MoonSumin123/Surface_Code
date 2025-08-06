import random
import pymatching
import sys
import matplotlib.pyplot as plt

class surface:
    def __init__(self):
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
    def __init__(self, 
                 x_boundary: bool = False, 
                 z_boundary: bool = False, 
                 t_x_boundary: bool = False, 
                 t_z_boundary: bool = False):

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

def error_making(list_of_data_qubit, probability):
    """
    Error make on data qubits of surface model with some probability

    Input:
        list_of_data_qubit: A list of data qubits of surface model
        probability: A int of physical error rate (Pauli error)

    Output:
        A list of data qubits that changed by some error rate
    """
    p = probability
    for data_qubit in list_of_data_qubit:
        data_qubit.value = random.choices([0,1,2,3], weights=[1-3*p, p, p, p], k=1)[0]

    return list_of_data_qubit

def measure_stabilizer(list_of_stabilizer, is_X=True):
    """
    Measure stabilizer's eigenvalue

    Input:
        list_of_stabilizer: A list of stabilizers of surface model (X or Z)

    Output:
        A list of stabilizers that measure their eigenvalue (check parity of neighbor data qubits error)
    """
    for stabilizer in list_of_stabilizer:
        stabilizer.value = 0  # 초기화 중요
        for data_qubit in stabilizer.data_qubits:
            if is_X:
                if data_qubit.value == 2 or data_qubit.value == 3:
                    stabilizer.value ^= 1
            else:
                if data_qubit.value == 1 or data_qubit.value == 3:
                    stabilizer.value ^= 1
    return list_of_stabilizer

def MWPM(list_of_stabilizer, is_X=True):
    """
    Minimum Weight Perfect Matching using pymatching

    Input:
        list_of_stabilizer: A list of stabilizers of surface model (X or Z)
        is_X: A bool: True if the stabilizer is of type X, else False (Z)

    Output:
        A n by 2 np.array containing n rows(edges) to be corrected, where each row is a tuple of two nodes(stabilizers).
    """
    pm = pymatching.Matching()
    
    visit = []
    syndrom = []
    for node in list_of_stabilizer:
        syndrom.append(node.value)
        visit.append(node)
        if node.is_boundary:
            pm.add_boundary_edge(node.index)
        for data_qubit in node.data_qubits:
            if is_X:
                neigbor_stabilizers = data_qubit.X_stabilizers
            else:
                neigbor_stabilizers = data_qubit.Z_stabilizers

            for node2 in neigbor_stabilizers:
                if node2 in visit:
                    continue
                pm.add_edge(node.index, node2.index)

    edge_info = pm.decode_to_edges_array(syndrom)
    return edge_info

def correction(list_of_stabilizer, edge_info, is_X=True):
    """
    Correct data qubits on surface by using the info of MWPM

    Input:
        list_of_stabilizer: A list of stabilizers of surface model (X or Z)
        edge_info: A n by 2 np.array that returned by the function 'MWPM'
        is_X: A bool: True if the stabilizer is of type X, else False (Z)

    Output:
        None
    """
    def flip(i, s):
        if s:
            mapping = {0:1, 1:0, 2:3, 3:2}
        else:
            mapping = {0:2, 2:0, 1:3, 3:1}
        return mapping[i]

    for edge in edge_info:
        index_node1, index_node2 = edge
        node1 = list_of_stabilizer[index_node1]
        if index_node2 == -1:
            for dq in node1.data_qubits:
                if is_X and dq.is_z_boundary:
                    dq.value = flip(dq.value, False)
                    break
                if not is_X and dq.is_x_boundary:
                    dq.value = flip(dq.value, True)
                    break

        else:
            node2 = list_of_stabilizer[index_node2]
            for data_qubit in node1.data_qubits:
                if is_X:
                    neigbor_stabilizers = data_qubit.X_stabilizers
                else:
                    neigbor_stabilizers = data_qubit.Z_stabilizers

                for neighbor_stabilizer in neigbor_stabilizers:
                    if neighbor_stabilizer == node2:
                        data_qubit.value = flip(data_qubit.value, not is_X)
                        
    return

def is_valid_correction(list_of_stabilizer):
    """
    Error correction is verified to be complete by confirming that all stabilizers are zero

    Input:
        list_of_stabilizer: A list of stabilizers of surface model (X or Z)

    Output:
        A bool: if True there is no non-zero, else False
    """
    for stabilizer in list_of_stabilizer:
        if stabilizer.value != 0:
            return False
        
    return True

def check_error_chain(list_of_data_qubit):
    """
    The error chain is identified by checking the boundary data qubits

    Input:
        list_of_data_qubit: A list of data qubits of surface model
        is_X: A bool: True if the stabilizer is of type X, else False (Z)

    Output:
        An int: 1 if there is at least one error chain, else 0
    """
    num_x_logical_error = 0
    num_z_logical_error = 0

    for dq in list_of_data_qubit:
        if dq.test_x_b:
            if dq.value in (1, 3):
                num_x_logical_error += 1
        if dq.test_z_b:
            if dq.value in (2, 3):
                num_z_logical_error += 1

    num_x_logical_error %= 2
    num_z_logical_error %= 2
    return 1 if num_x_logical_error or num_z_logical_error else 0

def one_cycle_simulation(surface, prob):
    """
    Simulate surface model one cycle (Collection of library functions)

    Input:
        surface: A class of model we want to simulate
        prob: A integer number the rate of error occurence

    Output:
        An int: 1 if there is at least one error chain, else 0
    """
    surface.data_qubits = error_making(surface.data_qubits, prob)
    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)
    edge_info_X = MWPM(surface.X_stabilizers, is_X=True)
    edge_info_Z = MWPM(surface.Z_stabilizers, is_X=False)

    # print("===== Data Qubits =====")
    # for i, dq in enumerate(surface.data_qubits):
    #     print(
    #         f"Qubit[{i}] val={dq.value}"
    #         f" | is_x_boundary={dq.is_x_boundary}"
    #         f" | is_z_boundary={dq.is_z_boundary}"
    #         f" | X_stabs={[x.index for x in dq.X_stabilizers]}"
    #         f" | Z_stabs={[z.index for z in dq.Z_stabilizers]}"
    #     )

    # print("edge_info_Z:", edge_info_Z)
    # print("edge_info_X:", edge_info_X)

    correction(surface.X_stabilizers, edge_info_X, is_X=True)
    correction(surface.Z_stabilizers, edge_info_Z, is_X=False)
    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)

    # print("===== Z Stabilizers =====")
    # for i, stab in enumerate(surface.Z_stabilizers):
    #     print(f"Z[{i}] val={stab.value} | boundary={stab.is_boundary} | data={[q.value for q in stab.data_qubits]}")

    # print("===== X Stabilizers =====")
    # for i, stab in enumerate(surface.X_stabilizers):
    #     print(f"X[{i}] val={stab.value} | boundary={stab.is_boundary} | data={[q.value for q in stab.data_qubits]}")

    # print("===== Data Qubits =====")
    # for i, dq in enumerate(surface.data_qubits):
    #     print(
    #         f"Qubit[{i}] val={dq.value}"
    #         f" | is_x_boundary={dq.is_x_boundary}"
    #         f" | is_z_boundary={dq.is_z_boundary}"
    #         f" | X_stabs={[x.index for x in dq.X_stabilizers]}"
    #         f" | Z_stabs={[z.index for z in dq.Z_stabilizers]}"
    #     )

    if not is_valid_correction(surface.X_stabilizers) or not is_valid_correction(surface.Z_stabilizers):
        print("Invalid Error Correction Occurred")
        raise RuntimeError("Invalid error correction at p={}".format(prob))

    logical_error = check_error_chain(surface.data_qubits)
    if logical_error:
        return 1
    return 0

def total_sumulation(list_of_distance, list_of_probability, iteration=1000):
    """
    Simulate surface model

    Input:
        list_of_distance: A list of distace to test
        list_of_probability: A list of probability to test
        iteration: A integer number for loop

    Output:
        None (print chart)
    """
    result_total = []   # [[d1, [[p1, rate11], [p2,rate12], ...]],
                        #  [d2, [[p1, rate21], [p2,rate22], ...]], ... ]
    for distance in list_of_distance:
        target_surface = surface(distance)
        result_distance = []
        for probability in list_of_probability:
            logical_error = 0
            for i in range(iteration):
                logical_error += one_cycle_simulation(target_surface, probability)
            rate_logical_error = logical_error / iteration
            result_distance.append([probability, rate_logical_error])
        result_total.append([distance, result_distance])

    make_chart(result_total)
    return

def make_chart(list_of_result):
    """
    Visulaizing chart

    Input:
        list_of_result: A list of result

    Output:
        None (print chart)
    """
    plt.figure(figsize=(8, 6))
    for distance_data in list_of_result:
        distance = distance_data[0]
        prob_rate_pairs = distance_data[1]
        probs = [pair[0] for pair in prob_rate_pairs]
        rates = [pair[1] for pair in prob_rate_pairs]
        plt.plot(probs, rates, marker='o', label=f'd = {distance}')

    plt.xlabel("Physical error rate (p)")
    plt.ylabel("Logical error rate")
    plt.title("Logical Error Rate vs Physical Error Rate")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return

