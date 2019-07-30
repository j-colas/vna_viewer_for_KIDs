#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import QWidget, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure

    
class MplWidget(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        
        #self.canvas.figure.set_visible(False)
        self.canvas.figure.set_facecolor("none")
        self.canvas.setStyleSheet("background-color:transparent;")
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        
        
        self.setLayout(vertical_layout)
        