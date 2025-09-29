import random
import pymatching
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
from my_class_ver1 import *

# --- BP-OSD helper (add to simulation.py) ---
try:
    from ldpc.bposd_decoder import BpOsdDecoder  # ldpc>=2.x
except Exception:
    from ldpc import BpOsdDecoder                # 일부 배포본 호환


def visualize_data_qubits_by_layer_subplot(surf):
    '''
    각 층별 data qubit의 Error 발생 여부 확인(시각화) 함수
    '''
    qubits = surf.qubits
    size = len(qubits)
    n_layers = size
    n_cols = math.ceil(math.sqrt(n_layers))
    n_rows = math.ceil(n_layers / n_cols)

    color_map = {0: 'lightgray', 1: 'red', 2: 'blue', 3: 'purple'}
    state = {0: 'None', 1: 'X', 2: 'Z', 3: 'Y'}
    legend_handles = [
        mpatches.Patch(color=color_map[v], label=f'Error={state[v]}')
        for v in sorted(color_map)
    ]

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3*n_cols, 3*n_rows),
                             squeeze=False)
    for layer in range(n_layers):
        row, col = divmod(layer, n_cols)
        ax = axes[row][col]

        xs, ys, vals = [], [], []
        for j in range(size):
            for k in range(size):
                obj = qubits[layer][j][k]
                if isinstance(obj, data_qubit):
                    xs.append(k)
                    ys.append(j)
                    vals.append(obj.value)

        ax.scatter(xs, ys,
                   c=[color_map[v] for v in vals],
                   s=80, marker='o')
        ax.set_title(f'Z = {layer}')
        # 눈금 켜기
        ax.set_xticks(range(size))
        ax.set_yticks(range(size))
        ax.set_xlim(-0.5, size-0.5)
        # 만약 (0,0)을 상단 왼쪽에 놓고 싶다면 아래 줄의 주석을 해제하세요
        # ax.set_ylim(size-0.5, -0.5)
        ax.set_aspect('equal')

    # 빈 칸은 숨기기
    total = n_rows * n_cols
    for empty in range(n_layers, total):
        row, col = divmod(empty, n_cols)
        axes[row][col].axis('off')

    fig.legend(handles=legend_handles,
               loc='lower center',
               ncol=len(color_map),
               title='data_qubit.value')
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()


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
        target_surface = build_rotated_surface_code(distance)
        result_distance = []
        for probability in list_of_probability:
            logical_error = 0
            p_eff = min(2*probability, 0.499)           # BP-OSD 확률 설정
            for i in range(iteration):
                logical_error += one_cycle_simulation(target_surface, probability, p_eff)
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


def one_cycle_simulation(surface, prob, p_eff):
    """
    Simulate surface model one cycle (Collection of library functions)

    Input:
        surface: A class of model we want to simulate
        prob: A integer number the rate of error occurence
        p_eff: A float number the probability of 실제 오류 (for BP-OSD의 신념 설정을 해)

    Output:
        An int: 1 if there is at least one error chain, else 0
    """
    surface.data_qubits = error_making(surface.data_qubits, prob)


    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)


    bits_X, cols_X, map_X = BPOSD_decode_bits(surface.X_stabilizers,
                                                    error_rate=p_eff, max_iter=50,
                                                    osd_order=2, osd_method='osd_cs')
    bits_Z, cols_Z, map_Z = BPOSD_decode_bits(surface.Z_stabilizers,
                                            error_rate=p_eff, max_iter=50,
                                            osd_order=2, osd_method='osd_cs')


    correction_bits(bits_X, cols_X, map_X, is_X=True)
    correction_bits(bits_Z, cols_Z, map_Z, is_X=False)


    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)


    if not is_valid_correction(surface.X_stabilizers) or not is_valid_correction(surface.Z_stabilizers):
        print("Invalid Error Correction Occurred")
        raise RuntimeError("Invalid error correction at p={}".format(prob))


    logical_error = check_error_chain(surface.data_qubits)
    if logical_error:
        return 1
    return 0


def error_making(list_of_data_qubit, probability, seed=None):
    """
    Error make on data qubits of surface model with some probability,
    with optional fixed random seed for reproducibility.

    Input:
        list_of_data_qubit: A list of data_qubit instances
        probability:        A float of physical error rate (Pauli error)
        seed:               An int or None. If int, RNG is seeded for reproducibility.

    Output:
        The same list_of_data_qubit, with each .value set according to the error model.
    """
    # 독립 RNG 인스턴스 생성 (seed=None 이면 system time 기반)
    rng = random.Random(seed)
    p = probability

    for data_qubit in list_of_data_qubit:
        # 0: no error, 1/2/3: Pauli errors X/Y/Z
        data_qubit.value = rng.choices([0, 1, 2, 3],
                                       weights=[1 - 3*p, p, p, p],
                                       k=1)[0]
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


def check_error_chain(list_of_data_qubit, surf = None):
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



# MWPM decoder & correction (for 2D surface code)
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
        if s:   # X flip
            mapping = {0:1, 1:0, 2:3, 3:2}
        else:   # Z flip
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


# BP-OSD decoder & correction (for 3D surface code)
def _build_pcm_and_cols(list_of_stabilizer):
    """
    안정자-데이터큐빗 인접성으로 패리티체크행렬 H 생성
    반환: H (m x k), cols (길이 k 의 '실제 dq.idx' 목록), idx_to_qubit (dq.idx -> 객체)
    """
    m = len(list_of_stabilizer)
    used = sorted({dq.idx for stab in list_of_stabilizer for dq in stab.data_qubits})
    col_index = {idx: i for i, idx in enumerate(used)}
    H = np.zeros((m, len(used)), dtype=np.uint8)        # parity matrix(행: dq, 열: stab)
    for stab in list_of_stabilizer:
        r = stab.index
        for dq in stab.data_qubits:
            H[r, col_index[dq.idx]] ^= 1
    idx_to_qubit = {dq.idx: dq for stab in list_of_stabilizer for dq in stab.data_qubits}
    return H, used, idx_to_qubit


def BPOSD_decode_bits(list_of_stabilizer, error_rate=0.1,
                      max_iter=50, osd_order=2, osd_method='osd_cs'):
    """
    syndrome -> 뒤집을 비트열(데이터큐빗 열들) 반환
    """
    H, cols, idx_to_qubit = _build_pcm_and_cols(list_of_stabilizer)     # parity matrix, 큐빗 정보, 큐빗index-object 매핑
    syndrome = np.array([stab.value for stab in list_of_stabilizer], dtype=np.uint8)
    dec = BpOsdDecoder(H,
                       error_rate=error_rate,
                       max_iter=max_iter,
                       bp_method='product_sum',
                       schedule='serial',
                       osd_method=osd_method,
                       osd_order=osd_order)
    bits = dec.decode(syndrome).astype(np.uint8)    # 1 이면 플립, 2 이면 그대로
    return bits, cols, idx_to_qubit                 # 플립 정보, 큐빗 정보, 큐빗index-object 매핑


def correction_bits(bits, cols, idx_to_qubit, is_X=True):
    """
    BP-OSD가 준 비트열에 따라 데이터 큐빗에 X/Z 플립 적용
    is_X=True: X 안정자 synd. => Z 플립 필요
    is_X=False: Z 안정자 synd. => X 플립 필요
    """
    def flip(i, do_X):
        if do_X:     # X flip
            return {0:1, 1:0, 2:3, 3:2}[i]
        else:        # Z flip
            return {0:2, 2:0, 1:3, 3:1}[i]

    do_X = (not is_X)  # X-stab에는 Z플립, Z-stab에는 X플립
    for j, b in enumerate(bits):
        if b:
            dq_idx = cols[j]
            dq = idx_to_qubit.get(dq_idx)
            if dq is not None:
                dq.value = flip(dq.value, do_X)


