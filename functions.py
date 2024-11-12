import numpy as np

mu = 0.3
_E = 2e11
t = 0.001

def read_ansys_nodes_and_elems():
    nodes = []
    elems = []
    
    with open( "nodes.txt", 'r') as file:
        next(file)
        
        for line in file:

            n, x, y, z = line.replace(",", ".").split()

            nodes.append([float(x) * 1000, -float(z) * 1000])
    
    with open("elems.txt", 'r') as file:
        next(file)
        
        for line in file:

            line = line.replace(",", ".").split()

            nodes_indexes = [int(ind) - 1 for ind in line[2:len(line)]]

            if len(nodes_indexes) == 3:
                elems.append(nodes_indexes)
            else:
                elems.append([nodes_indexes[0], nodes_indexes[1], nodes_indexes[2]])
                elems.append([nodes_indexes[3], nodes_indexes[2], nodes_indexes[0]])

    indexes = np.array([[3, 7, 8, 9, 10, 6], [25, 29, 28, 27, 26, 11]])
    lines = np.array([nodes[3], nodes[7], nodes[8], nodes[9], nodes[10], nodes[6]])

    return np.array(nodes), np.array(elems), lines, indexes

def gauss_solve(A, b):
    n = len(A)

    for i in range(n):
        max_idx = i

        for j in range(i + 1, n):
            if abs(A[j, i]) > abs(A[max_idx, i]):
                max_idx = j

        A[[i, max_idx]] = A[[max_idx, i]]
        b[[i, max_idx]] = b[[max_idx, i]]

        for j in range(i + 1, n):
            ratio = A[j, i] / A[i, i]
            A[j, i:] -= ratio * A[i, i:]
            b[j] -= ratio * b[i]

    X = np.zeros(n)

    for i in range(n - 1, -1, -1):
        X[i] = (b[i] - np.dot(A[i, i + 1:], X[i + 1:])) / A[i, i]

    return X

def get_E_matrix():
    return _E / (1 - mu * mu) * np.array([
            [1, mu, 0],
            [mu, 1, 0],
            [0, 0, (1 - mu)/2]
        ])

def calc_local_stiffness_matrices(elems, nodes):
    matrices = []
    Bm = []

    for elem in elems:
        k = np.zeros((6,6))
        nodes_elem = nodes[list(elem)]
        x, y = nodes_elem[:, 0], nodes_elem[:, 1]
        A = np.abs(0.5 * np.linalg.det([
            [1, x[0], y[0]],
            [1, x[1], y[1]],
            [1, x[2], y[2]]
        ]))

        B = 0.5 * A * np.array([
                [y[1] - y[2], 0, y[2] - y[0], 0,y[0] - y[1], 0],
                [0, x[2] - x[1], 0, x[0] - x[2], 0, x[1] - x[0]],
                [x[2] - x[1], y[1] - y[2], x[0] - x[2], y[2] - y[0], x[1] - x[0], y[0] - y[1]]
        ])

        E = get_E_matrix()

        k = t * A * np.dot(np.dot(B.T,E),B)

        matrices.append(k)
        Bm.append(B)

    return np.array(matrices), Bm

def calc_global_stiffness_matrix(elems, nodes):
    matrix = np.zeros((len(nodes) * 2, len(nodes) * 2))
    matrices, Bm = calc_local_stiffness_matrices(elems, nodes)

    for k in range(len(elems)):
        for i in range(6):
            for j in range(6):
                sm_x = 1 if j % 2 != 0 else 0
                sm_y = 1 if i % 2 != 0 else 0

                matrix[elems[k][i // 2] * 2 + sm_y][elems[k][j // 2] * 2 + sm_x] += matrices[k][i][j]

    return matrix, Bm

def make_F_vector(F):
    return np.array([[f[0] * np.cos(f[1]), f[0] * np.sin(f[1])] for f in F]).flatten()

def crossout_hinges_indexes(glob_matrix, H, F):
    for i in range(len(H)):
        if H[i] == 1:
            for j in range(len(glob_matrix)):
                glob_matrix[j][i] = 0 if i != j else 1
                glob_matrix[i][j] = 0 if i != j else 1

            F[i] = 0

    return glob_matrix

def calc_displacement(elems, nodes, F, H):
    glob_matrix, Bm = calc_global_stiffness_matrix(elems, nodes)

    Fm = make_F_vector(F)

    glob_matrix = crossout_hinges_indexes(glob_matrix, H, Fm)

    displacement = gauss_solve(glob_matrix, Fm).reshape(len(nodes), 2)

    Eps, Sig, Epsilons_eq, Sigmas_eq  = calc_deformation(elems, displacement, Bm)

    return displacement, Eps, Sig, Epsilons_eq, Sigmas_eq 

def calc_deformation(elems, d, Bm):
    Em = get_E_matrix()
    Bm = np.array(Bm)

    Epsilons = []
    Sigmas = []
    Epsilons_eq = []
    Sigmas_eq = []

    for i in range(len(elems)):
        Eps = Bm[i] @ d[elems[i]].flatten()
        Sigma = Eps @ Em

        e_x, e_y, e_xy = Eps
        s_x, s_y, s_xy = Sigma

        Eps_eq = (np.sqrt(2) / 3) * np.sqrt(((e_x - e_y) ** 2) + (e_y ** 2) + (e_x ** 2) + (1.5 * (e_xy ** 2)))
        Sig_eq = (1 / np.sqrt(2)) * np.sqrt(((s_x - s_y)) ** 2 + (s_y ** 2) + (s_x ** 2) + (6 * (s_xy ** 2)))

        Epsilons_eq.append(Eps_eq)
        Sigmas_eq.append(Sig_eq)
        Epsilons.append(Eps)
        Sigmas.append(Sigma)

    Sigmas = np.array(Sigmas)
    Epsilons = np.array(Epsilons)
    Epsilons_eq = np.array(Epsilons_eq)
    Sigmas_eq = np.array(Sigmas_eq)

    return Epsilons, Sigmas, Epsilons_eq, Sigmas_eq

def make_disposed_force(q, nodes):
    Fs = np.zeros(len(nodes))

    for i in range(0, len(nodes) - 1):
        l = np.abs(nodes[i][1] - nodes[i + 1][1])

        F = l * q / 2

        Fs[i] += F
        Fs[i + 1] += F

    return Fs

def make_disposed_F_vector(nodes, force_value, lines, indexes):
    F = make_disposed_force(force_value, lines)

    Fs = np.zeros((len(nodes), 2))
    
    i = 0

    for ind in indexes[0]:
        Fs[ind] = [F[i], np.pi]
        i += 1
    
    i = 0

    for ind in indexes[1]:
        Fs[ind] = [F[i], 0]
        i += 1
    
    return Fs

def make_hinges_supports(nodes, fixated):
    H = np.zeros(len(nodes) * 2)

    for fix in fixated:
        H[fix * 2] = 1
        H[fix * 2 + 1] = 1
    
    return H

def equation_of_line(start_point, end_point):
    x1, y1 = start_point
    x2, y2 = end_point

    if x2 - x1 != 0:
        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1

        return slope, intercept
    else:
        return float('inf'), x1

def intersection_point(line1, line2):
    m1, b1 = line1
    m2, b2 = line2

    if m1 == m2:
        return None

    if m1 == float('inf'):
        x = b1
        y = m2 * x + b2
    elif m2 == float('inf'):
        x = b2
        y = m1 * x + b1
    else:
        x = (b2 - b1) / (m1 - m2)
        y = m1 * x + b1

    return x, y

def split_area(horizontal_lines, vertical_lines, nodes_coords, ind):
    lines = []
    nodes = []

    for hline in horizontal_lines:
        line = []

        line1 = equation_of_line(hline[0], hline[1])

        for vline in vertical_lines:
            line2 = equation_of_line(vline[0], vline[1])

            x, y = intersection_point(line1, line2)
            
            y1 = np.abs(vline[0][1])
            y2 = np.abs(vline[1][1])

            y1 = min(y1, y2)
            y2 = max(y2, y1)

            if y >= vline[1][1] and y <= vline[0][1] and x >= hline[0][0] and x <= hline[1][0]:
                nodes_coords[(x, y)] = ind
                ind += 1

                line.append([x, y])
                nodes.append([x, y])
        
        lines.append(line)
    
    return np.array(lines), np.array(nodes), nodes_coords, ind

def make_elems(lines_from, lines_to, line_elem_from, line_elem_to, lines, nodes_coords):
    elems = []

    for i in range(lines_from, lines_to - 1):
        for j in range(line_elem_from, line_elem_to - 1):
            x1, y1 = lines[i][j]
            x2, y2 = lines[i][j + 1]
            x3, y3 = lines[i + 1][j + 1]
            x4, y4 = lines[i + 1][j]
            
            ind1 = nodes_coords[(x1, y1)]
            ind2 = nodes_coords[(x2, y2)]
            ind3 = nodes_coords[(x3, y3)]
            ind4 = nodes_coords[(x4, y4)]

            elems.append([ind1, ind2, ind3])
            elems.append([ind1, ind4, ind3])
    
    return np.array(elems)

def generate_extra_vertical_lines(x_from, x_to, y1, y2, div_size):
    div_step = (x_to - x_from) * div_size

    extra_lines_y = []

    for x in np.arange(x_from + div_step, x_to, div_step):
        extra_lines_y.append([[x, y1], [x, y2]])

    return extra_lines_y

def generate_extra_horizontal_lines(y_from, y_to, x1, x2, div_size):
    div_step = -np.abs((np.abs(y_from) - np.abs(y_to))) * div_size

    extra_lines_x = []

    for y in np.arange(y_from + div_step, y_to, div_step):

        extra_lines_x.append([[x1, y], [x2, y]])
    
    return extra_lines_x

def split_first_part(c, b, f, div_x_size, div_y_size, ind, nodes_coords):
    hlines = np.array([
        [[0, 0], [f, 0]],
        *generate_extra_horizontal_lines(0, -c, 0, f, div_y_size),
        [[0, -c], [f, -c]],
        *generate_extra_horizontal_lines(-c, -b, 0, f, div_y_size),
        [[0, -b], [f, -b]]
    ])

    vlines = np.array([
        [[0, 0], [0, -b]],
        *generate_extra_vertical_lines(0, f, 0, -b, div_x_size),
        [[f, 0], [f, -b]],
    ])
    
    return split_area(hlines, vlines, nodes_coords, ind)

def split_second_part(f, i, c, b, g, l, div_x_size, div_y_size, ind, nodes_coords):
    hlines = np.array([
        [[f + i, 0], [g, 0]],
        *generate_extra_horizontal_lines(0, -l, f + i, g, div_y_size),
        [[f + i, -c], [g, -l]],
    ])

    vlines = np.array([
        [[f + i, 0], [f + i, -c]],
        *generate_extra_vertical_lines(f + i, g, 0, -b, div_x_size),
         [[g, 0], [g, -b]]
    ])
    
    return split_area(hlines, vlines, nodes_coords, ind)

def split_third_part(c, l, b, j, f, i, g, div_x_size, div_y_size, ind, nodes_coords):
    hlines = np.array([
        [[f + i + i * 0.5, -c - c * 0.5], [g + j * 0.5, -l - c * 0.5]],
        *generate_extra_horizontal_lines(-c - c * 0.5, -b, f + i + i * 0.5, g + j * 0.5, div_y_size),
        [[f + i + i * 0.5, -b], [g - j * 0.5, -b]],
    ])

    vlines = np.array([
        [[f + i + i * 0.5, -c], [f + i + i * 0.5, -b]],
        *generate_extra_vertical_lines(f + i + i * 0.5, g - j * 0.5, -c, -b, div_x_size),
         [[g - j * 0.5, -l], [g - j * 0.5, -b]]
    ])
    
    return split_area(hlines, vlines, nodes_coords, ind)

def split_fourth_part(g, a, l, b, j, div_x_size, div_y_size, ind, nodes_coords):
    hlines = np.array([
        [[g + j, 0], [a, 0]],
        *generate_extra_horizontal_lines(0, -l, g + j, a, div_y_size),
        [[g + j, -l], [a, -l]],
        *generate_extra_horizontal_lines(-l, -b, g + j, a, div_y_size),
        [[g + j, -b], [a, -b]]
    ])

    vlines = np.array([
        [[g + j, 0], [g + j, -b]],
        *generate_extra_vertical_lines(g + j, a, 0, -b, div_x_size),
        [[a, 0], [a, -b]],
    ])
    
    return split_area(hlines, vlines, nodes_coords, ind)

def split_grid(a, b, f, g, c, l, i, j, div_vertical_size = 0.5, div_horizontal_size = 0.5):
    nodes_coords = {}
    ind = 0
    elems = []

    lines1, nodes1, nodes_coords, ind = split_first_part(c, b, f, div_vertical_size, div_horizontal_size, ind, nodes_coords)
    lines2, nodes2, nodes_coords, ind = split_second_part(f, i, c, b, g, l, div_vertical_size * 0.5, div_horizontal_size, ind, nodes_coords)
    lines3, nodes3, nodes_coords, ind = split_third_part(c, l, b, j,  f, i, g, div_vertical_size * 0.5, div_horizontal_size, ind, nodes_coords)
    lines4, nodes4, nodes_coords, ind = split_fourth_part(g, a, l, b, j, div_vertical_size , div_horizontal_size, ind, nodes_coords)
    
    lines23 = [lines2[len(lines2) - 1], lines3[0]]

    c_ind = len(lines3)
    l_ind = len(lines3)

    lines_i = [lines1[-c_ind:][:,-1], lines3[:,0]]
    lines_j = [lines3[:,-1], lines4[-l_ind:][:,0]]

    elems1 = make_elems(0, len(lines1), 0, len(lines1[0]), lines1, nodes_coords)
    elems2 = make_elems(0, len(lines2), 0, len(lines2[0]), lines2, nodes_coords)
    elems3 = make_elems(0, len(lines3), 0, len(lines3[0]), lines3, nodes_coords)
    elems4 = make_elems(0, len(lines4), 0, len(lines4[0]), lines4, nodes_coords)
    elems23 = make_elems(0, len(lines23), 0, len(lines23[0]), lines23, nodes_coords)
    elems_i = make_elems(0, len(lines_i), 0, len(lines_i[0]), lines_i, nodes_coords)
    elems_j = make_elems(0, len(lines_j), 0, len(lines_j[0]), lines_j, nodes_coords)

    x1, y1 = lines1[-c_ind:][:,-1][0]
    x2, y2 = lines2[-1][0]
    x3, y3 = lines3[0][0]

    elems.append([nodes_coords[(x1, y1)], nodes_coords[(x2, y2)], nodes_coords[(x3, y3)]])

    lines3[:,-1], lines4[-l_ind:][:,0]

    x1, y1 = lines4[-l_ind:][:,0][0]
    x2, y2 = lines2[-1][-1]
    x3, y3 = lines3[0][-1]

    elems.append([nodes_coords[(x1, y1)], nodes_coords[(x2, y2)], nodes_coords[(x3, y3)]])

    elems = np.concatenate((elems, elems1, elems2, elems3, elems4, elems23, elems_i, elems_j))
    nodes = np.concatenate((nodes1, nodes2, nodes3, nodes4))

    first_indexes = [nodes_coords[(x, y)] for x, y in lines1[:][:,0]]
    last_indexes = [nodes_coords[(x, y)] for x, y in lines4[:][:,-1]]

    hinges_indexes = first_indexes
    force_indexes = [first_indexes, last_indexes]

    indexes = [first_indexes, last_indexes]
    lines = lines1[:][:,0]

    return elems, nodes, indexes, lines