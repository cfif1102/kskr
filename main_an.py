import functions as funcs
import graphics as graphs
import numpy as np

q = 10e3

nodes, elems, lines, indexes = funcs.read_ansys_nodes_and_elems()

F = funcs.make_disposed_F_vector(nodes, q, lines, indexes)
H = funcs.make_hinges_supports(nodes, indexes[0])

#graphs.plot_fig(elems, nodes, F, H)

d, eps, sig, epsilons_eq, sigmas_eq = funcs.calc_displacement(elems, nodes, F, H)

new_nodes = nodes + d
print(np.sort(epsilons_eq)[::-1][:15])
#graphs.plot_fig(elems, new_nodes, F, H)

graphs.show_equal_deformations(new_nodes, elems, epsilons_eq, sigmas_eq)
