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
        self.qubits = []

    def add_data(self, q):
        self.data_qubits.append(q)

    def add_x(self, x):
        self.X_stabilizers.append(x)

    def add_z(self, z):
        self.Z_stabilizers.append(z)

class data_qubit:
    def __init__(self, idx: int, x_boundary: bool = False, z_boundary: bool = False, t_x_boundary: bool = False, t_z_boundary: bool = False):
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

        self.test_x_b = False        # test for X logical Error (이름 잘못 정해놨음. X 오류 확인하려는 목적이 맞음)
        self.test_z_b = False        # test for Z logical Error (이름 잘못 정해놨음. Z 오류 확인하려는 목적이 맞음)

        self.idx = idx

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
    """
    매트릭스를 통해 3D cubic lattice 만드는 함수
    """
    surf = surface()
    size = 2 * d - 1

    qubits = [[[None for _ in range(size)] for _ in range(size)] for _ in range(size)]

    qubit_idx = 0
    for i in range(size):                                   # 층
        if (i%2) == 0:
            for j in range(size):                           # 행
                for k in range(size):                       # 열
                    if (j%2) == 0 and (k%2)==0:             # Data quvit
                        q = data_qubit(qubit_idx)
                        qubits[i][j][k] = q
                        surf.add_data(q)
                        qubit_idx += 1
                        if j == 0 and k == 0:               # test for X logical error
                            q.test_x_b = True
                        if i == 0:                          # test for Z logical error
                            q.test_z_b = True
                    elif (j%2)==1 and (k%2)==1:             # None
                        continue
                    else:                                   # Z stabilizer
                        stab = stabilizer()
                        stab.set_ind(len(surf.Z_stabilizers))
                        surf.add_z(stab)
                        qubits[i][j][k] = stab
        else:
            for j in range(size):                           # 행
                for k in range(size):                       # 열
                    if (j%2) == 0 and (k%2)==0:             # X stabilizer
                        stab = stabilizer()
                        stab.set_ind(len(surf.X_stabilizers))
                        surf.add_x(stab)
                        qubits[i][j][k] = stab
                    elif (j%2)==1 and (k%2)==1:             # Z stabilizer
                        stab = stabilizer()
                        stab.set_ind(len(surf.Z_stabilizers))
                        surf.add_z(stab)
                        qubits[i][j][k] = stab
                    else:                                   # Data qubit
                        q = data_qubit(qubit_idx)
                        qubits[i][j][k] = q
                        surf.add_data(q)
                        qubit_idx += 1

    neighbor_deltas = [
        (+1,  0,  0),  # 아래층
        (-1,  0,  0),  # 윗층
        ( 0, +1,  0),  # 행 아래
        ( 0, -1,  0),  # 행 위
        ( 0,  0, +1),  # 열 오른쪽
        ( 0,  0, -1),  # 열 왼쪽
    ]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                q = qubits[i][j][k]
                # data_qubit만 처리
                if not isinstance(q, data_qubit):
                    continue

                # 6방향 이웃 탐색
                for di, dj, dk in neighbor_deltas:
                    ni, nj, nk = i + di, j + dj, k + dk

                    # 배열 경계 체크
                    if not (0 <= ni < size and 0 <= nj < size and 0 <= nk < size):
                        continue

                    neighbor = qubits[ni][nj][nk]
                    # stabilizer가 있을 때만 연결
                    if isinstance(neighbor, stabilizer):
                        # stabilizer 쪽에도 이 data_qubit 추가
                        neighbor.add_data(q)

                        # X, Z 중 어느 리스트에 들어있는지 보고 data_qubit에도 추가
                        if neighbor in surf.X_stabilizers:
                            q.add_x(neighbor)
                        elif neighbor in surf.Z_stabilizers:
                            q.add_z(neighbor)

                # data qubit의 boundary 지정하는 코드 (현재 simulation code에서는 사용하지 않음. 만약 MWPM을 사용한다면 다시 사용해야 핳 듯)
                # if len(q.X_stabilizers) < 2:
                #     q.is_z_boundary = True
                # if len(q.Z_stabilizers) < 4:
                #     q.is_x_boundary = True

    # stabilizer의 boundary 지정하는 코드 (현재 simulation code에서는 사용하지 않음. 만약 MWPM을 사용한다면 다시 사용해야 핳 듯)
    # for x_stab in surf.X_stabilizers:
    #     x_stab.is_boundary = any(q.is_z_boundary for q in x_stab.data_qubits)
    # for z_stab in surf.Z_stabilizers:
    #     z_stab.is_boundary = any(q.is_x_boundary for q in z_stab.data_qubits)

    surf.qubits = qubits
    return surf

