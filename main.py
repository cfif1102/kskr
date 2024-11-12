import graphics as graph
import functions as funcs   
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QVBoxLayout,  QMainWindow, QVBoxLayout, QPushButton, QGroupBox, QHBoxLayout, QSlider, QVBoxLayout, QMessageBox,QToolBar
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QDesktopServices
import numpy as np
import numpy as np

a = 15
b = 5
f = 3
g = 12
c = 2
l = 2
i = 0.5
j = 0.5

class PlotWidget(QWidget):
    def __init__(self):
        super().__init__()

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.figure = Figure(figsize=(15, 15))
        self.canvas = FigureCanvas(self.figure)

        layout.addWidget(self.canvas)

    def add_subplot(self, *args, **kwargs):
        return self.figure.add_subplot(*args, **kwargs)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
    
        self.grid_vals = [0.2, 0.25, 0.5, 1.0]

        self.d = None
        self.eps = None
        self.sig = None
        self.eps_eq = None
        self.sig_eq = None
        self.is_splitted = False
        self.elems = None
        self.new_nodes = None
        self.is_splitted = False

        self.setWindowTitle("Анализ элемента трубопровода")

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        layout = QHBoxLayout(self.central_widget)

        self.plot_widget = PlotWidget()
        self.plot_layout = QVBoxLayout()

        self.plot_layout.addWidget(self.plot_widget)

        layout.addLayout(self.plot_layout)

        row_layout = QVBoxLayout()

        #-----------------------------------
        division_group_box = QGroupBox("Разбиение фигуры")
        col_division_box = QVBoxLayout(division_group_box)
        division_group_box.setLayout(col_division_box) 

        row_elems_layout = QHBoxLayout()

        self.divsion_step_x_text = 'Шаг сетки по оси X: '
        divsion_step_x_label = QLabel(self.divsion_step_x_text + "1", self)
        self.divsion_step_x_label = divsion_step_x_label
        divsion_step_x = QSlider(Qt.Horizontal, self)

        self.divsion_step_x = divsion_step_x
        self.divsion_step_x.setFocusPolicy(Qt.NoFocus)
        self.divsion_step_x.setRange(0, len(self.grid_vals) - 1)
        self.divsion_step_x.setSingleStep(1)
        self.divsion_step_x.valueChanged[int].connect(self.onSliderChanged1)
        self.divsion_step_x.setValue(len(self.grid_vals) - 1)

        row_elems_layout.addWidget(divsion_step_x_label)
        row_elems_layout.addWidget(divsion_step_x)

        col_division_box.addLayout(row_elems_layout)
        #------------------------

        row_elems_layout = QHBoxLayout()

        self.divsion_step_y_text = 'Шаг сетки по оси Y: '
        divsion_step_y_label = QLabel(self.divsion_step_y_text + "1", self)
        self.divsion_step_y_label = divsion_step_y_label
        divsion_step_y = QSlider(Qt.Horizontal, self)

        self.divsion_step_y = divsion_step_y
        self.divsion_step_y.setFocusPolicy(Qt.NoFocus)
        self.divsion_step_y.setRange(0, len(self.grid_vals) - 1)
        self.divsion_step_y.setSingleStep(1)
        self.divsion_step_y.valueChanged[int].connect(self.onSliderChanged2)
        self.divsion_step_y.setValue(len(self.grid_vals) - 1)

        row_elems_layout.addWidget(divsion_step_y_label)
        row_elems_layout.addWidget(divsion_step_y)

        col_division_box.addLayout(row_elems_layout)

        #-----------------------

        row_elems_layout = QHBoxLayout()

        self.split_figure_btn = QPushButton("Построить сетку", self)
        self.split_figure_btn.clicked.connect(self.split_figure_action)

        row_elems_layout.addWidget(self.split_figure_btn)

        col_division_box.addLayout(row_elems_layout)

        #----------------------

        row_layout.addWidget(division_group_box)

        division_group_box = QGroupBox("Параметры")
        col_division_box = QVBoxLayout(division_group_box)
        division_group_box.setLayout(col_division_box) 

        #----------------------

        row_elems_layout = QHBoxLayout()

        self.q_label = QLabel('Распределенная нагрузка q (В ньютонах):', self)
        self.q_value = QLineEdit(self)

        row_elems_layout.addWidget(self.q_label)
        row_elems_layout.addWidget(self.q_value)

        col_division_box.addLayout(row_elems_layout)

        #----------

        #----------

        row_elems_layout = QHBoxLayout()

        self.puasson_ration_label = QLabel('Коэффициент Пуассона (0.1 - 0.5):', self)
        self.puasson_ration_value = QLineEdit(self)

        row_elems_layout.addWidget(self.puasson_ration_label)
        row_elems_layout.addWidget(self.puasson_ration_value)

        col_division_box.addLayout(row_elems_layout)

        #----------

        row_elems_layout = QHBoxLayout()

        self.thickness_label = QLabel('Толщина:', self)
        self.thickness = QLineEdit(self)

        row_elems_layout.addWidget(self.thickness_label)
        row_elems_layout.addWidget(self.thickness)

        col_division_box.addLayout(row_elems_layout)

        #----------

        row_elems_layout = QHBoxLayout()

        self.make_calculations_btn = QPushButton("Провести расчеты", self)
        self.make_calculations_btn.clicked.connect(self.make_calculations)

        row_elems_layout.addWidget(self.make_calculations_btn)

        col_division_box.addLayout(row_elems_layout)

        #----------

        row_layout.addWidget(division_group_box)

        division_group_box = QGroupBox("Построение графиков")
        col_division_box = QVBoxLayout(division_group_box)
        division_group_box.setLayout(col_division_box) 

        #----- Графики

        row_elems_layout = QHBoxLayout()

        self.stress_btn = QPushButton("Напряжение", self)
        self.stress_btn.clicked.connect(self.plot_sigmas)

        row_elems_layout.addWidget(self.stress_btn)

        col_division_box.addLayout(row_elems_layout)

        #--------

        row_elems_layout = QHBoxLayout()

        self.strain_btn = QPushButton("Деформация", self)
        self.strain_btn.clicked.connect(self.plot_epsilons)

        row_elems_layout.addWidget(self.strain_btn)

        col_division_box.addLayout(row_elems_layout)

        #--------

        row_elems_layout = QHBoxLayout()

        self.eq_stress = QPushButton("Эквивалентное напряжение и деформация", self)
        self.eq_stress.clicked.connect(self.plot_eq_sigmas)

        row_elems_layout.addWidget(self.eq_stress)

        col_division_box.addLayout(row_elems_layout)

        #--------

        row_elems_layout = QHBoxLayout()

        self.clear_field = QPushButton("Очистить поле", self)
        self.clear_field.clicked.connect(self.clean_window)

        row_elems_layout.addWidget(self.clear_field)

        col_division_box.addLayout(row_elems_layout)

        #--------

        row_layout.addWidget(division_group_box)

        layout.addLayout(row_layout)

        elems, nodes, indexes, lines = funcs.split_grid(a, b, f, g, c, l, i, j, 1, 1)

        self.c_elems = elems
        self.c_nodes = nodes

        graph.plot_fig(elems, nodes, plot_widget = self.plot_widget)
    
    def onSliderChanged1(self, value):
        selected_value = self.grid_vals[value]

        self.div_x = float(self.grid_vals[value])

        self.divsion_step_x_label.setText(self.divsion_step_x_text + str(selected_value))
    
    def onSliderChanged2(self, value):
        selected_value = self.grid_vals[value]

        self.div_y = float(self.grid_vals[value])

        self.divsion_step_y_label.setText(self.divsion_step_y_text + str(selected_value))
    
    def plot_eq_sigmas(self):
        if self.can_be_ploted():
            graph.show_equal_deformations(self.new_nodes, self.elems, self.eps_eq, self.sig_eq, self.plot_widget)
        else:
            QMessageBox.warning(self, 'Ошибка', f'Проведите расчеты для построения графиков!')

    def plot_epsilons(self):
        if self.can_be_ploted():
            self.is_eps_selected = True
            graph.plot_deformations(self.new_nodes, self.elems, self.eps, "Epsilon", self.plot_widget)
        else:
            QMessageBox.warning(self, 'Ошибка', f'Проведите расчеты для построения графиков!')
    
    def plot_sigmas(self):
        if self.can_be_ploted():
            graph.plot_deformations(self.new_nodes, self.elems, self.sig, "Sigma", self.plot_widget)
        else:
            QMessageBox.warning(self, 'Ошибка', f'Проведите расчеты для построения графиков!')

    def set_null_vals(self):
        self.d = None
        self.eps = None
        self.sig = None
        self.eps_eq = None
        self.sig_eq = None
        self.is_splitted = False
        self.elems = None
        self.new_nodes = None
        self.is_splitted = False

    def make_calculations(self):
        if self.is_splitted == False:
            QMessageBox.warning(self, 'Ошибка', f'Перед проведением расчетов постройте сетку!')
            return

        try:
            q_value = float(self.q_value.text())
            puasson_ration_value = float(self.puasson_ration_value.text())
            thickness_value = float(self.thickness.text())

            if q_value < 0 or puasson_ration_value < 0 or thickness_value < 0:
                QMessageBox.warning(self, 'Ошибка', f'Параметры не могут быть отрицательными!')
                return

            if puasson_ration_value < 0.1 or puasson_ration_value > 0.5:
                QMessageBox.warning(self, 'Ошибка', f'Коэффициент Пуассона должен быть в пределах от 0.1 до 0.5!')
                return

            funcs._mu = puasson_ration_value
            funcs.t = thickness_value


            elems, nodes, indexes, lines = funcs.split_grid(a, b, f, g, c, l, i, j, self.div_x, self.div_y)

            F = funcs.make_disposed_F_vector(nodes, q_value, lines, indexes)
            H = funcs.make_hinges_supports(nodes, indexes[0])

            displacement, Eps, Sig, Epsilons_eq, Sigmas_eq = funcs.calc_displacement(elems, nodes, F, H)

            self.d = displacement
            self.eps = Eps
            self.sig = Sig
            self.eps_eq = Epsilons_eq
            self.sig_eq = Sigmas_eq
            self.new_nodes = nodes + displacement
            self.elems = elems

            graph.plot_fig(self.elems, self.new_nodes, F, H, plot_widget = self.plot_widget)
        except:
            QMessageBox.warning(self, 'Ошибка', f'Проверьте введенные данные!')
            return

    def can_be_ploted(self):
        return self.d is not None
    
    def clean_window(self):
        self.set_null_vals()

        graph.plot_fig(self.c_elems, self.c_nodes, plot_widget = self.plot_widget)

    def split_figure_action(self):
        div_x = self.grid_vals[int(self.divsion_step_x.value())]
        div_y = self.grid_vals[int(self.divsion_step_y.value())]

        self.set_null_vals()
        
        self.is_splitted = True

        elems, nodes, indexes, lines = funcs.split_grid(a, b, f, g, c, l, i, j, div_x, div_y)

        graph.plot_fig(elems, nodes, plot_widget = self.plot_widget)

app = QApplication(sys.argv)

window = MainWindow()
window.show()

sys.exit(app.exec_())


'''
a = 15
b = 5
f = 3
g = 12
c = 2
l = 2
i = 0.5
j = 0.5

q = 10e3

elems, nodes, indexes, lines = funcs.split_grid(a, b, f, g, c, l, i, j, 0.2, 0.25)

F = funcs.make_disposed_F_vector(nodes, q, lines, indexes)
H = funcs.make_hinges_supports(nodes, indexes[0])


#graphs.plot_fig(elems, nodes, F, H)


d, eps, sig, epsilons_eq, sigmas_eq = funcs.calc_displacement(elems, nodes, F, H)

print(np.sort(epsilons_eq)[::-1][:15])

new_nodes = nodes + d

#graphs.plot_fig(elems, new_nodes, F, H)

graphs.plot_deformations(new_nodes, elems, eps, sig)


graphs.show_equal_deformations(new_nodes, elems, epsilons_eq, sigmas_eq)
'''