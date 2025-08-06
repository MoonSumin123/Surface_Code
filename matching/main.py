from test import *

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