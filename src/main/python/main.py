#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from fbs_runtime.application_context.PyQt5 import ApplicationContext
from PyQt5.QtWidgets import QMainWindow, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import numpy as np
import ntpath
import sys
import calibKiD as cK

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

class MatplotlibWidget(QMainWindow):
    
    def __init__(self):
        
        QMainWindow.__init__(self)

        loadUi(appctxt.get_resource("qt_designer.ui"),self)
        #loadUi('qt_designer.ui',self)
        self.setWindowTitle("VNA Scan Viewer For KIDs")
        self.setWindowIcon(QIcon(appctxt.get_resource("icon.png")))
        
        self.pushButton_generate_random_signal.setEnabled(False)
        self.fitButton.setEnabled(True)
        self.headerButton.setEnabled(False)
        self.browseI.setEnabled(False)
        self.browseQ.setEnabled(False)
        self.pushButton_generate_random_signal.clicked.connect(self.update_graph)
        self.fitButton.clicked.connect(self.add_fit)
        self.headerButton.clicked.connect(self.load_header)
        self.addToolBar(NavigationToolbar(self.MplWidget.canvas, self))
        
       
        self.browseCalib.clicked.connect(self.browse)
        self.browseI.clicked.connect(self.browse)
        self.browseQ.clicked.connect(self.browse)
        self.browseHead.clicked.connect(self.browse)
        
        self.clearButton.clicked.connect(self.clear_graph)
        
        #self.labelTone.editingFinished.connect(self.check_enable)
        #self.fnameCalib.textChanged.connect(self.check_enable)
        
        self.fcalib = None
        self.Icalib = None
        self.Qcalib = None
        self.tone = None
        self.run_number = None 
        
        self.watchdog = False
              

#    def check_enable(self):
#        
#        if self.fnameCalib.text() != "" and self.labelTone.text() != "":
#            self.tone = np.float32(self.labelTone.text())
#            self.enable_fit(self.fitButton)
#            
#    def enable_fit(self,Button):
#        Button.setEnabled(True)
        

    def load_header(self):
        params = cK.get_info(self.fnameHead.text(),self.fnameCalib.text())
        self.run_number = params[0]
        self.tone = params[6]
        self.labelRun.setText(self.run_number)
        
        str_tone = '%.0f'% self.tone
        self.labelTone.setText(str_tone)
        #self.check_enable(self)
        
    def browse(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","Text Files (*.txt);;Binary Files (*.bin);;All Files (*)", options=options)
        sender = self.sender()
        if fileName:
            print(fileName)
            if sender == self.browseCalib:
                self.fnameCalib.setText(fileName)
                self.pushButton_generate_random_signal.setEnabled(True)
            elif sender == self.browseI:
                self.fnameI.setText(fileName)
            elif sender == self.browseQ:
                self.fnameQ.setText(fileName)
            else :
                self.fnameHead.setText(fileName)
                self.headerButton.setEnabled(True)
                #self.check_enable(self)
                
    def clear_graph(self):
        self.MplWidget.canvas.axes.clear()
        self.MplWidget_2.canvas.axes.clear()
        self.MplWidget_3.canvas.axes.clear()
        self.MplWidget.canvas.draw()
        self.MplWidget_2.canvas.draw()
        self.MplWidget_3.canvas.draw()
            
    def update_graph(self):

        fname = self.fnameCalib.text()
        runName = path_leaf(fname)[4:11:]
        self.labelRun.setText(runName)
        
        data = np.genfromtxt(fname,comments="%")
        
        
        f = data[::,0]
        I = data[::,1]
        Q = data[::,2]
        
        self.fcalib = f
        self.Icalib = I
        self.Qcalib = Q
        
        normS21 = np.sqrt(I**2+Q**2)
        self.tone = f[np.argmin(normS21)]
        
        self.labelTone.setText('%.0f'% self.tone)
        
        #self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(I, Q, marker='o',color='k',alpha=.3,ls='',ms=5)

        #self.MplWidget.canvas.axes.legend(('IQ'),loc='upper right')
        self.MplWidget.canvas.axes.set_title('IQ Circle')
        self.MplWidget.canvas.axes.set_aspect('equal')
        self.MplWidget.canvas.draw()
        
        #self.MplWidget_2.canvas.axes.clear()
        self.MplWidget_2.canvas.axes.plot(f,normS21, marker='o',color='k',alpha=.3,ls='',ms=5)
        self.MplWidget_2.canvas.axes.set_title('Resonance')
        self.MplWidget_2.canvas.draw()
        
        #self.MplWidget_3.canvas.axes.clear()
        self.MplWidget_3.canvas.axes.plot(f,np.arctan2(I,Q), marker='o',color='k',alpha=.3,ls='',ms=5)
        self.MplWidget_3.canvas.axes.set_title('Phase')
        self.MplWidget_3.canvas.draw()
        
    def add_fit(self):
        if self.Icalib is None or self.tone is None:
            if not self.watchdog:
                self.update_graph()
                self.add_fit()
                self.watchdog = True
        else:      
            #print('plot fit')
            self.tone = np.float32(self.labelTone.text())
            self.labelTone.setText('%.0f'% self.tone)
            # f,I,Q,tone = fIQ
            params = cK.fit_S21_model([self.fcalib,self.Icalib,self.Qcalib,self.tone])
            [fr,Qc,Ql,phi0,tau,a,alpha,Pz,Rz,offset,xc,yc,r] = params
            Qi = 1./(1./Ql-1./Qc)
            
            fr_txt = '%03d'% fr
            self.labelFres.setText(fr_txt)
            Qc_txt = '%03d'% int(Qc/1000)
            self.labelQc.setText(Qc_txt)
            Ql_txt = '%03d'% int(Ql/1000)
            self.labelQl.setText(Ql_txt)
            Qi_txt = '%03d'% int(Qi/1000)
            self.labelQi.setText(Qi_txt)
            phi0_txt = '%05.3f'% phi0
            self.labelPhi0.setText(phi0_txt)
            tau_txt = '%05.3f'% tau
            self.labelTau.setText(tau_txt)
            a_txt = '%05.3f'% a
            self.labelA.setText(a_txt)
            alpha_txt = '%05.3f'% alpha
            self.labelAlpha.setText(alpha_txt)
            
            f_fit = np.linspace(self.fcalib.min(),self.fcalib.max(),2000)
            S21_fit = a*np.exp(1j*alpha)*cK.fun_model(f_fit,Qc,fr,Ql,phi0)
    
            i_res = np.argmin(np.abs(f_fit-fr))
            S21_fit_res = S21_fit[i_res]
            #idmin = np.argmin(np.abs(f_fit-fr))
            #X = np.complex(S21_fit[idmin].real,S21_fit[idmin].imag)
            #Pz = Pz*a*np.exp(1j*alpha)
    
            self.MplWidget.canvas.axes.plot(xc,yc, marker='o',mfc='deepskyblue',mec='k',color='r',alpha=.9,ls='-',ms=5,lw=2)
            self.MplWidget.canvas.axes.plot(np.real(S21_fit), np.imag(S21_fit), marker='',color='r',alpha=.7,ls='-',ms=1,lw=2)
            self.MplWidget.canvas.axes.plot(np.real(S21_fit_res), np.imag(S21_fit_res), marker='o',mfc='gold',mec='k',color='r',alpha=.9,ls='',ms=7,lw=2)        
            self.MplWidget.canvas.draw()
            
            self.MplWidget_2.canvas.axes.plot(f_fit, np.abs(S21_fit), marker='',color='r',alpha=.7,ls='-',ms=1,lw=2)
            self.MplWidget_2.canvas.axes.plot(f_fit[i_res],np.abs(S21_fit_res), marker='o',mfc='gold',mec='k',color='r',alpha=.9,ls='',ms=7,lw=2)
            self.MplWidget_2.canvas.draw()
            
            self.MplWidget_3.canvas.axes.plot(f_fit, np.arctan2(S21_fit.real,S21_fit.imag), marker='',color='r',alpha=.7,ls='-',ms=1,lw=2)
            self.MplWidget_3.canvas.axes.plot(f_fit[i_res],np.arctan2(S21_fit_res.real,S21_fit_res.imag), marker='o',mfc='gold',mec='k',color='r',alpha=.9,ls='',ms=7,lw=2)
            self.MplWidget_3.canvas.draw()

if __name__ == '__main__':
    appctxt = ApplicationContext()
    window = MatplotlibWidget()
    window.show()
    exit_code = appctxt.app.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)
