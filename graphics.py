import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap

def plot_grid(nodes, lines_vertical, lines_horizontal, a, b, f, g, c, l, i, j):
    plt.figure(figsize=(20, 10))

    plt.plot(nodes[:,0], nodes[:,1], 'ro', color="black")
    plt.fill([0, f, f, f + i, f + i, g, g, g + j, g + j, a, a, 0, 0], [0, 0, -c, -c, 0, 0, -l, -l, 0, 0, -b, -b, 0], color='gray', alpha=0.65)

    for i in range(len(nodes)):
        x, y = nodes[i]

        #plt.text(x + 0.1, y + 0.2, str(i), fontsize = 10, ha='left', va='center')

    for line in lines_vertical:
        x1, y1 = line[0]
        x2, y2 = line[1]

        plt.plot([x1, x2], [y1, -b], color='black')

    for line in lines_horizontal:
        x1, y1 = line[0]
        x2, y2 = line[1]

        plt.plot([x1, a], [y1, y2], color='black')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('off')
    plt.axis('equal')
    plt.margins(x=0.2, y=0.2)
    plt.grid(True)
    plt.show()

def plot_fig(elems, nodes, F = [], H = [], plot_widget = None):
    plot_widget.figure.clear()
    plot_widget.figure.tight_layout()

    ax = plot_widget.add_subplot(111)

    if len(H) > 0:
        indexes = np.where(H == 1)
        indexes = [index // 2 for index in indexes]
        nodes_locked = nodes[indexes][0]

        ax.plot(nodes_locked[:, 0], nodes_locked[:, 1], 'ro', color = 'green')

    for i in range(len(nodes)):
        x, y = nodes[i]

        ax.text(x + 0.1, y + 0.1, str(i), fontsize = 7, ha='left', va='center')

    for i in range(len(elems)):
        x = nodes[elems[i], 0]
        y = nodes[elems[i], 1]

        x = np.append(x, x[0])
        y = np.append(y, y[0])

        ax.fill(x, y, color = 'gray', edgecolor="black", alpha=0.85)

    if len(F) > 0:
        for i in range(len(F)):
            if F[i][0] != 0:
                x_end = 0.7 * np.cos(F[i][1])
                y_end = 0.7 * np.sin(F[i][1])

                x_start = nodes[i][0]
                y_start = nodes[i][1]

                ax.arrow(x_start, y_start, x_end, y_end, head_width=0.1, head_length=0.2, fc='red', ec='red')

    ax.axis('equal')
    ax.axis('off')
    ax.margins(x=0.05, y=0.05)
    ax.grid(False)

    plot_widget.canvas.draw()

def plot_deformations(nodes, elems, vals, title, plot_widget = None):
    plot_widget.figure.clear()
    plot_widget.figure.tight_layout()

    colors = ['#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000']
    cmap = LinearSegmentedColormap.from_list('custom', colors, N=256)

    Titles = ['_x', '_y', '_xy']

    axs = [plot_widget.add_subplot(1, 3, 1), plot_widget.add_subplot(1, 3, 2), plot_widget.add_subplot(1, 3, 3)]

    for i in range(3):
        eps = vals[:, i]

        triang = Triangulation(nodes[:, 0], nodes[:, 1], elems)

        norm = Normalize(vmin=eps.min(), vmax=eps.max())

        pcm = axs[i].tripcolor(triang, eps, cmap=cmap, norm=norm)

        axs[i].set_title(title + Titles[i])
        axs[i].axis('equal')
        axs[i].grid(True)
        axs[i].axis('off')

    plot_widget.figure.tight_layout()
    plot_widget.canvas.draw()

def show_plt(nodes, elems, values, ax, title, cmap):
    triang = Triangulation(nodes[:, 0], nodes[:, 1], elems)

    norm = Normalize(vmin=values.min(), vmax=values.max())

    pcm = ax.tripcolor(triang, values, cmap=cmap, norm=norm)
    plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8, aspect=20)

    ax.set_title(title)
    ax.margins(x=0.2, y=0.2 )
    ax.axis('equal')
    ax.axis('off')
    ax.grid(True)

def show_equal_deformations(nodes, elems, epsilons_eq, sigmas_eq, plot_widget = None):
    plot_widget.figure.clear()
    plot_widget.figure.tight_layout()

    colors = ['#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000']
    cmap = LinearSegmentedColormap.from_list('custom', colors, N=256)

    ax1, ax2 = [plot_widget.add_subplot(1, 2, 1), plot_widget.add_subplot(1, 2, 2)]

    show_plt(nodes, elems, epsilons_eq, ax1, 'Эквивалентная деформация', cmap)
    show_plt(nodes, elems, sigmas_eq, ax2, 'Эквивалентное напряжение', cmap)

    plot_widget.canvas.draw()