import random


class Qubit:
    name = "Qubit"
    
    def __init__(self):
        self.neighbors = []
        self.num_neighbor = 0
    
    def link(self, other):
        if not isinstance(other, Qubit):
            raise TypeError
        
        if other not in self.neighbors:
            self.neighbors.append(other)
            self.num_neighbor += 1
        if self not in other.neighbors:
            other.neighbors.append(self)
            other.num_neighbor += 1

    def unlink(self, other):
        if not isinstance(other, Qubit):
            raise TypeError
        
        try:
            self.neighbors.remove(other)
            self.num_neighbor-= 1
            other.neighbors.remove(self)
            other.num_neighbor -= 1
        except ValueError as e:
            print("There is no element in Qubit.neighbors list: ", e)

    def overlap(self, other):
        if not isinstance(other, Qubit):
            raise TypeError
        for o_neighbor in other.neighbors.copy():
            if o_neighbor in self.neighbors:
                other.unlink(o_neighbor)
            else:
                self.link(o_neighbor)
                other.unlink(o_neighbor)
        return other.num_neighbor == 0


class DataQubit(Qubit):
    name = "Data Qubit"
    error_rate = 0.3
    error_types = ["No Error", "X Error", "Z Error", "Y Error"]
    
    def __init__(self):
        super().__init__()
        self.error_type = "No Error"
        self.negative_sign = False

    def _toggle_stabilizers(self, neighbor_type):
        for stabilizer in self.neighbors:
            if isinstance(stabilizer, neighbor_type):
                stabilizer.toggle()

    def toggle_Xerror(self):
        if self.error_type == "No Error":
            self.error_type = "X Error"
        elif self.error_type == "X Error":
            self.error_type = "No Error"
        elif self.error_type == "Z Error":
            self.error_type = "Y Error"
            self.negative_sign = not self.negative_sign
        elif self.error_type == "Y Error":
            self.error_type = "Z Error"
            self.negative_sign = not self.negative_sign
        self._toggle_stabilizers(ZStabilizer)

    def toggle_Zerror(self):
        if self.error_type == "No Error":
            self.error_type = "Z Error"
        elif self.error_type == "X Error":
            self.error_type = "Y Error"
        elif self.error_type == "Z Error":
            self.error_type = "No Error"
        elif self.error_type == "Y Error":
            self.error_type = "X Error"
        self._toggle_stabilizers(XStabilizer)

    def toggle_Yerror(self):
        if self.error_type == "No Error":
            self.error_type = "Y Error"
        elif self.error_type == "X Error":
            self.error_type = "Z Error"
        elif self.error_type == "Z Error":
            self.error_type = "X Error"
            self.negative_sign = not self.negative_sign
        elif self.error_type == "Y Error":
            self.error_type = "No Error"
            self.negative_sign = not self.negative_sign
        self._toggle_stabilizers((ZStabilizer, XStabilizer))

    def next_cycle(self):
        r = random.random()
        if (r < DataQubit.error_rate/3):
            self.toggle_Xerror()
        elif (r < DataQubit.error_rate*2/3):
            self.toggle_Zerror()
        elif (r < DataQubit.error_rate):
            self.toggle_Yerror()
        else:
            pass
        return self.error_type

    def reset(self):
        self.error_type = "No Error"
        self.negative_sign = False

class UStabilizer(Qubit):
    name = "Universal Stabilizer Qubit"
    def __init__(self):
        super().__init__()
        self.syndrome = False
    def toggle(self):
        self.syndrome = not self.syndrome
    def reset(self):
        self.syndrome = False

class XStabilizer(UStabilizer):
    name = "X Stabilizer Qubit"
    def __init__(self):
        super().__init__()

class ZStabilizer(UStabilizer):
    name = "Z Stabilizer Qubit"
    def __init__(self):
        super().__init__()



class SurfaceGraph:
    def __init__(self):
        self.data_qubit = []
        self.X_stabilizer = []
        self.Z_stabilizer = []
    
    def make_graph(self):
        pass


    class UBoundaryFigure:
        def __init__(self, t: str, k: int):
            self.t = t
            self.k = k
            self.make_figure(t, k)
            self.overlapped_figures = [None for i in range(k)]
        
        def make_figure(self, t, k):
            if t == "X":
                self.center = ZStabilizer()
                self.stabilizers = [XStabilizer() for i in range(k)]
            elif t == "Z":
                self.center = XStabilizer()
                self.stabilizers = [ZStabilizer() for i in range(k)]
            else:
                raise ValueError
            self.data_qubits = [DataQubit() for i in range(k)]
            
            prev_stab_qubit = self.stabilizers[-1]
            for d_qubit, stab_qubit in zip(self.data_qubits, self.stabilizers):
                self.center.link(d_qubit)
                prev_stab_qubit.link(d_qubit)
                stab_qubit.link(d_qubit)
                prev_stab_qubit = stab_qubit

        def overlap_same_figure(self, other, index_s):
            if self.t != other.t or self.k != other.k:
                raise TypeError
            while index_s < 0:
                index_s += self.k
            while index_s >= self.k:
                index_s -= self.k
            index_o = index_s + other.k // 2
            if index_o >= other.k:  index_o -= other.k
            stab1 = self.stabilizers[index_s]
            stab2 = other.staabilizers[index_o]
            stab1.overlap(stab2)
            other.stabilizers[index_o] = stab1


    class ZBoundaryFigureTypeA(UBoundaryFigure):
        def __init__(self):
            super().__init__("Z", 7)

    class ZBoundaryFigureTypeB(UBoundaryFigure):
        def __init__(self):
            super().__init__("Z", 3)

    class XBoundaryFigureTypeC(UBoundaryFigure):
        def __init__(self):
            super().__init__("X", 4)
