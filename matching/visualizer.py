import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.colors import ListedColormap

from test import *

def visualize_surface_code_correction(surface, probability, show_steps=True):
    if show_steps:
        # Step 0: 초기 상태 (에러 없음)
        visualize_surface_state(surface, "Step 0: Initial State (No Errors)", step=0)
    
    # Step 1: 에러 발생
    surface.data_qubits = error_making(surface.data_qubits, probability)
    if show_steps:
        visualize_surface_state(surface, f"Step 1: After Physical Errors (p={probability})", step=1)
    
    # Step 2: Syndrome 측정
    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)
    if show_steps:
        visualize_surface_state(surface, "Step 2: After Syndrome Measurement", step=2)
    
    # Step 3: MWPM으로 에러 경로 찾기
    edge_info_X = MWPM(surface.X_stabilizers, is_X=True)
    edge_info_Z = MWPM(surface.Z_stabilizers, is_X=False)
    if show_steps:
        visualize_surface_with_matching(surface, edge_info_X, edge_info_Z, "Step 3: MWPM Correction Paths", step=3)
    
    # Step 4: 에러 정정 수행
    correction(surface.X_stabilizers, edge_info_X, is_X=True)
    correction(surface.Z_stabilizers, edge_info_Z, is_X=False)
    if show_steps:
        visualize_surface_state(surface, "Step 4: After Error Correction", step=4)
    
    # Step 5: 최종 검증
    surface.X_stabilizers = measure_stabilizer(surface.X_stabilizers, is_X=True)
    surface.Z_stabilizers = measure_stabilizer(surface.Z_stabilizers, is_X=False)
    logical_error = check_error_chain(surface.data_qubits)
    
    if show_steps:
        title = f"Step 5: Final State - {'LOGICAL ERROR' if logical_error else 'SUCCESS'}"
        visualize_surface_state(surface, title, step=5, final=True, logical_error=logical_error)
    
    return logical_error

def visualize_surface_state(surface, title, step=0, final=False, logical_error=False):
    """
    Surface Code의 현재 상태를 시각화
    """
    d = int(np.sqrt(len(surface.data_qubits)) + 1) // 2 + 1
    size = 2 * d - 1
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    data_positions = {}
    stabilizer_positions_X = {}
    stabilizer_positions_Z = {}
    
    # Data qubit 위치 계산
    data_idx = 0
    for i in range(size):
        for j in range(size):
            if (i + j) % 2 == 0:
                if data_idx < len(surface.data_qubits):
                    data_positions[(i, j)] = surface.data_qubits[data_idx]
                    data_idx += 1
    
    # Stabilizer 위치 계산
    x_stab_idx = 0
    z_stab_idx = 0
    for i in range(size):
        for j in range(size):
            if i % 2 == 0 and j % 2 == 1:  # X stabilizer
                if x_stab_idx < len(surface.X_stabilizers):
                    stabilizer_positions_X[(i, j)] = surface.X_stabilizers[x_stab_idx]
                    x_stab_idx += 1
            elif i % 2 == 1 and j % 2 == 0:  # Z stabilizer
                if z_stab_idx < len(surface.Z_stabilizers):
                    stabilizer_positions_Z[(i, j)] = surface.Z_stabilizers[z_stab_idx]
                    z_stab_idx += 1
    
    # Data qubit 그리기
    for (i, j), qubit in data_positions.items():
        color = get_qubit_color(qubit.value)
        circle = patches.Circle((j, size-1-i), 0.3, facecolor=color, edgecolor='black', linewidth=2)
        ax.add_patch(circle)
        
        error_text = get_error_text(qubit.value)
        if error_text:
            ax.text(j, size-1-i, error_text, ha='center', va='center', fontsize=10, fontweight='bold')
        
        if qubit.is_x_boundary or qubit.is_z_boundary:
            boundary_circle = patches.Circle((j, size-1-i), 0.35, fill=False, edgecolor='red', linewidth=3)
            ax.add_patch(boundary_circle)
    
    # X Stabilizer
    for (i, j), stab in stabilizer_positions_X.items():
        color = 'red' if stab.value == 1 else 'lightcoral'
        alpha = 1.0 if stab.value == 1 else 0.3
        rect = patches.Rectangle((j-0.2, size-1-i-0.2), 0.4, 0.4, 
                               facecolor=color, alpha=alpha, edgecolor='darkred', linewidth=2)
        ax.add_patch(rect)
        
        if stab.value == 1:
            ax.text(j, size-1-i, 'X', ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    
    # Z Stabilizer
    for (i, j), stab in stabilizer_positions_Z.items():
        color = 'blue' if stab.value == 1 else 'lightblue'
        alpha = 1.0 if stab.value == 1 else 0.3
        rect = patches.Rectangle((j-0.2, size-1-i-0.2), 0.4, 0.4, 
                               facecolor=color, alpha=alpha, edgecolor='darkblue', linewidth=2)
        ax.add_patch(rect)
        
        if stab.value == 1:
            ax.text(j, size-1-i, 'Z', ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    
    ax.set_xlim(-0.5, size-0.5)
    ax.set_ylim(-0.5, size-0.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    title_color = 'red' if (final and logical_error) else 'green' if final else 'black'
    ax.set_title(title, fontsize=14, fontweight='bold', color=title_color)

    add_legend(ax, step)
    
    plt.tight_layout()
    plt.show()

def visualize_surface_with_matching(surface, edge_info_X, edge_info_Z, title, step=3):
    """
    MWPM 결과와 함께 Surface Code 상태 시각화
    """
    visualize_surface_state(surface, title, step)
    
    fig = plt.gcf()
    ax = plt.gca()
    
    d = int(np.sqrt(len(surface.data_qubits)) + 1) // 2 + 1
    size = 2 * d - 1
    
    # X stabilizer matching 경로 (빨간 선)
    for edge in edge_info_X:
        if len(edge) == 2:
            idx1, idx2 = edge
            if idx2 != -1:  # boundary가 아닌 경우
                pos1 = get_stabilizer_position(surface.X_stabilizers[idx1], size, is_X=True)
                pos2 = get_stabilizer_position(surface.X_stabilizers[idx2], size, is_X=True)
                if pos1 and pos2:
                    ax.plot([pos1[1], pos2[1]], [pos1[0], pos2[0]], 'r-', linewidth=3, alpha=0.7)
    
    # Z stabilizer matching 경로 (파란 선)
    for edge in edge_info_Z:
        if len(edge) == 2:
            idx1, idx2 = edge
            if idx2 != -1:  # boundary가 아닌 경우
                pos1 = get_stabilizer_position(surface.Z_stabilizers[idx1], size, is_X=False)
                pos2 = get_stabilizer_position(surface.Z_stabilizers[idx2], size, is_X=False)
                if pos1 and pos2:
                    ax.plot([pos1[1], pos2[1]], [pos1[0], pos2[0]], 'b-', linewidth=3, alpha=0.7)
    
    plt.show()

def get_stabilizer_position(stabilizer, size, is_X=True):
    """
    Stabilizer의 grid 상에서의 위치를 찾는 함수
    """
    stab_idx = 0
    for i in range(size):
        for j in range(size):
            if is_X and i % 2 == 0 and j % 2 == 1:
                if stab_idx == stabilizer.index:
                    return (size-1-i, j)
                stab_idx += 1
            elif not is_X and i % 2 == 1 and j % 2 == 0:
                if stab_idx == stabilizer.index:
                    return (size-1-i, j)
                stab_idx += 1
    return None

def get_qubit_color(value):
    colors = {
        0: 'white',      # I (no error)
        1: 'lightblue',  # Z error
        2: 'lightcoral', # X error  
        3: 'plum'        # Y error (X+Z)
    }
    return colors.get(value, 'gray')

def get_error_text(value):
    texts = {
        0: '',
        1: 'Z',
        2: 'X', 
        3: 'Y'
    }
    return texts.get(value, '?')

def add_legend(ax, step):
    legend_elements = [
        patches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', label='No Error'),
        patches.Circle((0, 0), 0.1, facecolor='lightblue', edgecolor='black', label='Z Error'),
        patches.Circle((0, 0), 0.1, facecolor='lightcoral', edgecolor='black', label='X Error'),
        patches.Circle((0, 0), 0.1, facecolor='plum', edgecolor='black', label='Y Error'),
    ]
    
    if step >= 2:  # Syndrome 측정 이후
        legend_elements.extend([
            patches.Rectangle((0, 0), 0.1, 0.1, facecolor='red', label='X Syndrome'),
            patches.Rectangle((0, 0), 0.1, 0.1, facecolor='blue', label='Z Syndrome'),
        ])
    
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))

def demo_visualization():
    surf = build_rotated_surface_code(3)
    logical_error = visualize_surface_code_correction(surf, 0.1, show_steps=True)
    
    print(f"Logical Error Occurred: {logical_error}")

surf = build_rotated_surface_code(3)
logical_error = visualize_surface_code_correction(surf, 0.1, show_steps=True)