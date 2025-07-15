import random
import pymatching
import sys
import matplotlib.pyplot as plt



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

# 이 코드가 문제. 가장 쉽게는 추가적으로 확인하는 boundary qubit list를 만들면 되긴 한데..