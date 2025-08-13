from __future__ import print_function, division

"""
    Written by Weverson R. Gomes, 2015.
    Copyright (c) 2015, Weverson R. Gomes.

    This file is part of Sasha Software.

    Sasha Software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sasha Software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sasha Software.  If not, see <http://www.gnu.org/licenses/>.
"""

# This is only needed for Python v2 but is harmless for Python v3.
import sip
sip.setapi('QString', 2)
sip.setapi('QVariant', 2)

#from matplotlib.backends import qt4_compat
#use_pyside = qt4_compat.QT_API == qt4_compat.QT_API_PYSIDE

#if use_pyside:
#    from PySide.QtCore import QtGui, QtCore
#    from PySide.QtGui import QThread, SIGNAL
#else:
#    from PyQt4 import QtGui, QtCore
#    from PyQt4.QtCore import QThread, SIGNAL

from PyQt4.QtGui import QApplication,QPixmap,QSplashScreen,QProgressBar, QAction
from PyQt4.QtCore import Qt
import time
import sys
import os
import parameters

# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    exedir= os.path.dirname(sys.executable)
elif __file__:
    exedir = os.path.dirname(os.path.abspath(__file__))
#exedir = os.path.dirname(os.path.realpath(__file__))

class MyProgressBar(QProgressBar):
    def __init__(self, parent = None):
        QProgressBar.__init__(self, parent)
        self.setStyleSheet(parameters.stylesheetqt("styleProgressBar"))

    def setValue(self, value):
        QProgressBar.setValue(self, value)

        if value == self.maximum():
            self.setStyleSheet(parameters.stylesheetqt("styleData"))

# SPLASH SCREEN
app = QApplication(sys.argv)
splash_pix = QPixmap(exedir+'/images/aboutlogo.png')
splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)

# adding progress bar
progressBar = MyProgressBar(splash)
progressBar.setGeometry(splash.width()/11, 8*splash.height()/8.8,8*splash.width()/10, splash.height()/16)

# SHOWING SPLASH + PROGRESS BAR
splash.setMask(splash_pix.mask())
splash.show()

i0 = 0
progressBar.setValue(i0)
from PyQt4.QtGui import  QFileDialog,QFileDialog,QMessageBox,QMainWindow,QPrinter,QAbstractPrintDialog,QMovie,QDialog,QDialogButtonBox,QPainter,QPrinter,QPrintDialog,QPrintPreviewDialog
i0 += 1
progressBar.setValue(i0)
from PyQt4.QtCore import QObject,QThread,SIGNAL,pyqtSignal,QCoreApplication,Qt,QRect,QFileInfo,QFile,QTextStream,QSettings,QMetaObject
i0 += 1
progressBar.setValue(i0)
from logging import Handler,getLogger,Formatter,DEBUG
i0 += 1
progressBar.setValue(i0)
import os
i0 += 1
progressBar.setValue(i0)
from HF import hartree_fock
i0 += 1
progressBar.setValue(i0)
import design
i0 += 1
progressBar.setValue(i0)
from get_basis import inp
i0 += 1
progressBar.setValue(i0)
from collections import Counter
i0 += 1
progressBar.setValue(i0)
from itertools import repeat
i0 += 1
progressBar.setValue(i0)


PACKAGENAME = 'SASHA PROGRAM'



from matplotlib import use
use('Qt4Agg')
i0 += 1
progressBar.setValue(i0)
from matplotlib.pyplot import close,rcParams,figure,legend,ylabel,xlabel,tick_params,xlim
i0 += 1
progressBar.setValue(i0)

# PRINT OUPUT LOG IN QT
class QtHandler(Handler):
    def __init__(self):
        Handler.__init__(self)
    def emit(self, record):
        record = self.format(record)
        if record: XStream.stdout().write('%s\n'%record)
        # originally: XStream.stdout().write("{}\n".format(record))

logger = getLogger(__name__)
handler = QtHandler()
handler.setFormatter(Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(DEBUG)


class XStream(QObject):
    _stdout = None
    _stderr = None
    messageWritten = pyqtSignal(str)
    def flush( self ):
        pass
    def fileno( self ):
        return -1
    def write( self, msg ):
        if ( not self.signalsBlocked() ):
            self.messageWritten.emit(msg)
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
        return XStream._stderr

class ConsoleWindowLogHandler(Handler):
    def __init__(self, sigEmitter):
        super(ConsoleWindowLogHandler, self).__init__()
        self.sigEmitter = sigEmitter

    def emit(self, logRecord):
        message = str(logRecord.getMessage())
        self.sigEmitter.emit(SIGNAL("logMsg(QString)"), message)


class getThread(QThread):
    '''
    Make a new thread instance for calculation
    '''

    def __init__(self,N,R,Z,A,D,CartAng,Center,n_electrons,n,lenL,locL,lenAtomN,atoms,mass,MOstart,MOend,method,basis,charge,units,inputName,fnCube,dictcheckBox,plot_type):
        QThread.__init__(self)
        self.N = N
        self.R = R
        self.Z = Z
        self.A = A
        self.D = D
        self.CartAng = CartAng
        self.Center = Center
        self.n_electrons = n_electrons
        self.n = n
        self.lenL = lenL
        self.locL = locL
        self.lenAtomN = lenAtomN
        self.atoms = atoms
        self.mass = mass
        self.MOstart = MOstart
        self.MOend = MOend
        self.method = method
        self.basis = basis
        self.charge = charge
        self.units = units
        self.inputName = inputName
        self.fnCube = fnCube
        self.dictcheckBox = dictcheckBox
        self.plot_type = plot_type

    def __del__(self):
        self.quit()

    def stop(self):
        while self.isRunning():
            self.terminate()
            QCoreApplication.processEvents()

    def calculation(self):

        global e, count, countlist, EN, ENT, ENlist, P, C, errorcalculation
        errorcalculation = False
        count = False
        e, count, countlist,  EN, ENT, ENlist, P, C, errorcalculation = hartree_fock(self.N,self.R,self.Z,self.A,self.D,self.CartAng,self.Center,self.n_electrons,self.n,self.lenL,self.locL,self.lenAtomN,self.atoms,self.mass,self.MOstart,self.MOend,self.method,self.basis,self.charge,self.units,self.inputName,self.fnCube,self.dictcheckBox,self.plot_type,logger)

    def run(self):
        '''
        Go to main SCF function
        '''
        run = self.calculation()

class sashaprogram(QMainWindow, design.Ui_MainWindow):
    """
    How the basic structure of PyQt GUI code looks and behaves like is
    explained in this tutorial
    http://nikolak.com/pyqt-qt-designer-getting-started/
    """

    def __init__(self, fileName=None):
        super(self.__class__, self).__init__()
        self.setupUi(self)


        self.aboutDialog = None
        self.locationDialog = None
        self.fnCube = False
        self.plot_type = False

        # LINK MENU ACTION WITH ACTIONS
        self.connect(self.action_Open, SIGNAL("triggered()"), self.fileOpen)
        self.connect(self.action_Save, SIGNAL("triggered()"), self.fileSave)
        self.connect(self.action_SaveAs, SIGNAL("triggered()"), self.fileSaveAs)
        self.connect(self.action_SaveOutputAs, SIGNAL("triggered()"), self.fileSaveOutputAs)
        self.connect(self.action_Print, SIGNAL("triggered()"), self.filePrint)
        self.connect(self.action_PrintPreview, SIGNAL("triggered()"), self.filePrintPreview)
        self.connect(self.action_About, SIGNAL("triggered()"), self.about)
        self.connect(self.actionExit, SIGNAL("triggered()"), self.maybeSave)

        self.runButton.clicked.connect(self.checkrun)
        self.orbitalButton.clicked.connect(self.plotorbital)
        self.convergenceButton.clicked.connect(self.plotconvergence)
        self.MODensityButton.clicked.connect(self.MODensityCall)
        self.printoutputButton.clicked.connect(self.fileSaveOutputAs)
        

        if fileName is None:
            fileName = ":/input.inp"
            self.fileName = ":/input.inp"

        self.setCurrentFileName()
        #self.updateRecentFileActions()

        # TEST IF INPUT WAS MODIFIED
        self.inputTextEdit.document().modificationChanged.connect(self.action_SaveAs.setEnabled)
        self.inputTextEdit.document().modificationChanged.connect(self.setWindowModified)
        self.setWindowModified(self.inputTextEdit.document().isModified())
        self.inputTextEdit.document().modificationChanged.connect(self.action_Save.setEnabled)
        self.action_Save.setEnabled(self.inputTextEdit.document().isModified())


        XStream.stdout().messageWritten.connect( self.outputTextEdit.insertPlainText )
        XStream.stderr().messageWritten.connect( self.outputTextEdit.insertPlainText )
        

    def checkrun(self):

        self.outputTextEdit.clear()

        self.basis = self.comboBox_2.currentText()

        self.method = self.comboBox.currentText()

        self.charge = self.chargeSpinBox.value()

        self.units = self.unitsComboBox.currentText()

        if self.maybeRun() and self.testinp():
            self.tabs.setCurrentIndex(3)

            self.dictcheckBox = [self.checkBox.isChecked(),self.checkBox_2.isChecked(),self.checkBox_3.isChecked(),self.checkBox_4.isChecked(),self.checkBox_5.isChecked(),self.checkBox_6.isChecked()]
            #dictcheckBox = {self.checkBox:False,self.checkBox_2:False,self.checkBox_3:False,self.checkBox_4:False,self.checkBox_5:False,self.checkBox_6:False}
            #for ck,test in dictcheckBox.items():
            #    if ck.isChecked():
            #        dictcheckBox[ck] = True

            #self.dictcheckBox = dictcheckBox.values()

            self.N,self.R,self.Z,self.A,self.D,self.CartAng,self.Center,self.n_electrons,self.n,self.lenL,self.locL,self.lenAtomN,self.atoms,self.mass = inp(self.basis,self.units,self.fileName)
            self.MOstart = self.MOend = False

            if self.checkBox_6.isChecked():

                if self.MOquestion():

                    self.fnCube = QFileDialog.getSaveFileName(self, "Save Cube File as...",exedir+"/"+self.fileName.split('/')[-1].split('.')[0],
                             "Cube File (*.cube);;All Files (*)",QFileDialog.DontUseNativeDialog)

                    if not self.fnCube:
                        self.checkBox_6.setChecked(False)
                        self.dictcheckBox = [self.checkBox.isChecked(),self.checkBox_2.isChecked(),self.checkBox_3.isChecked(),self.checkBox_4.isChecked(),self.checkBox_5.isChecked(),self.checkBox_6.isChecked()]
                    else:
                        lfn = self.fnCube.lower()
                        if not lfn.endswith(('.cube')):
                            # The default.
                            self.fnCube += '.cube'

                else:
                    self.checkBox_6.setChecked(False)
                    self.dictcheckBox = [self.checkBox.isChecked(),self.checkBox_2.isChecked(),self.checkBox_3.isChecked(),self.checkBox_4.isChecked(),self.checkBox_5.isChecked(),self.checkBox_6.isChecked()]

            #print(self.dictcheckBox)
            self.initialize()

        elif not self.testinp():
            self.tabs.setCurrentIndex(2)
            self.inputerror()

        else:
            pass


    def initialize(self):

        self.runHF = getThread(self.N,self.R,self.Z,self.A,self.D,self.CartAng,self.Center,self.n_electrons,self.n,self.lenL,self.locL,self.lenAtomN,self.atoms,self.mass,self.MOstart,self.MOend,self.method,self.basis,self.charge,self.units,self.fileName,self.fnCube,self.dictcheckBox,self.plot_type)

        # SIGNALS WHEN PUSH RUN AND KILL BUTTUN, LINKED WITH runHF THREAD
        self.connect(self.runHF, SIGNAL("started()"), self.sigstart)
        self.connect(self.runHF, SIGNAL("finished()"), self.sigdone)
        self.connect(self.runHF, SIGNAL("terminated()"), self.sigdoneterminated)

        #dummyEmitter = QObject()
        #self.connect(dummyEmitter, SIGNAL("logMsg(QString)"),self.outputTextEdit.appendPlainText)
        #self.connect(dummyEmitter, SIGNAL("logMsg(QString)"),self.outputTextEdit.append)
        #consoleHandler = ConsoleWindowLogHandler(dummyEmitter)
        #logger.setLevel(logging.INFO)
        #logger.addHandler(consoleHandler)

        self.runHF.start()

        #self.killButton.clicked.connect(self.runHF.stop)
        self.connect(self.killButton, SIGNAL("clicked()"), self.runHF.stop)

        import subprocess

        #proc = subprocess.Popen(['python', 'test.py'], shell=False)

    def MOquestion(self):
        #dialogMO = QDialog()
        dialogMO = design.MODialog(self.N,self.n_electrons)
        #dialogMO.ui.setupUi(dialogMO,self.N)
        #dialogMO.setAttribute(Qt.WA_DeleteOnClose)
        if dialogMO.exec_():
            self.MOstart = dialogMO.comboBox_start.currentText()
            self.MOend = dialogMO.comboBox_end.currentText()
            self.plot_type = dialogMO.comboBox_MO_DENSITY.currentText()
            return True
        else:
            return False

        #self.MO = PopupCube()
        #self.MO.setGeometry(QRect(100, 100, 400, 200))
        #self.MO.show()


    #@pyqtSlot(str)
    #def on_logBuffer_bufferMessage(self, message):
    #    self.outputTextEdit.append(message)

    def add_post(self,post_text):
        self.outputTextEdit.append(post_text)

    def sigstart(self):
        """
        Show the message.
        Disable Run button, enable the Stop one
        """
        self.killButton.setEnabled(True)
        self.runButton.setEnabled(False)
        self.orbitalButton.setEnabled(False)
        self.convergenceButton.setEnabled(False)
        self.MODensityButton.setEnabled(False)
        self.printoutputButton.setEnabled(False)
        self.statusBar().showMessage('Hartre Fock Calculation Was Started')

    def sigdone(self):
        """
        Show the message.
        Disable Stop button, enable the Start one
        """
        self.killButton.setEnabled(False)
        self.runButton.setEnabled(True)
        self.action_SaveOutputAs.setEnabled(True)
        self.tabs.setCurrentIndex(3)
        if self.dictcheckBox[4]:
            self.f = open(self.fileName.split('/')[-1].split('.')[0]+'_'+self.basis+'_'+'charge_'+str(self.charge)+'.out','w')
            self.fileSaveOutput()
        if count:
            import numpy as np
            self.ENplot =  np.asarray(ENlist)
            self.countplot = np.asarray(countlist)
            self.statusBar().showMessage('Hartre Fock Calculation Finished Successfully')
            self.orbitalButton.setEnabled(True)
            self.convergenceButton.setEnabled(True)
            self.MODensityButton.setEnabled(True)
            self.printoutputButton.setEnabled(True)

        elif errorcalculation:
            self.statusBar().showMessage('Hartre Fock Calculation Finished With Error')
            self.orbitalButton.setEnabled(False)
            self.convergenceButton.setEnabled(False)
            self.MODensityButton.setEnabled(False)
            self.printoutputButton.setEnabled(False)

        else:
            self.statusBar().showMessage('Hartre Fock Calculation Stopped')
            self.orbitalButton.setEnabled(False)
            self.convergenceButton.setEnabled(False)
            self.MODensityButton.setEnabled(False)
            self.printoutputButton.setEnabled(False)

        self.checkBox.setChecked(False)
        self.checkBox_2.setChecked(False)
        self.checkBox_3.setChecked(False)
        self.checkBox_4.setChecked(False)
        self.checkBox_5.setChecked(False)
        self.checkBox_6.setChecked(False)

    def sigdoneterminated(self):
        """
        Show the message.
        Disable Stop button, enable the Start one
        """
        self.killButton.setEnabled(False)
        self.runButton.setEnabled(True)
        self.orbitalButton.setEnabled(False)
        self.MODensityButton.setEnabled(False)
        self.printoutputButton.setEnabled(False)
        #self.convergenceButton.setEnabled(False)
        self.statusBar().showMessage('Stopping...')

        self.checkBox.setChecked(False)
        self.checkBox_2.setChecked(False)
        self.checkBox_3.setChecked(False)
        self.checkBox_4.setChecked(False)
        self.checkBox_5.setChecked(False)
        self.checkBox_6.setChecked(False)

    def setCurrentFileName(self, fileName=''):
        self.fileName = fileName
        self.inputTextEdit.document().setModified(False)

        if not fileName:
            shownName = 'untitled.inp'
        else:
            shownName = QFileInfo(fileName).fileName()
       
        settings = QSettings(PACKAGENAME, 'Recent Files')
        files = settings.value('recentFileList', [])



        if fileName:
            files.insert(0, fileName)


        settings.setValue('recentFileList', files)

        for widget in QApplication.topLevelWidgets():
            if isinstance(widget, sashaprogram):
                widget.updateRecentFileActions()

        self.setWindowTitle(self.tr("%s[*] - %s" % (shownName, PACKAGENAME)))
        self.setWindowModified(False)

    def fileNew(self):
        if self.maybeSave():
            self.inputTextEdit.clear()
            self.outputTextEdit.clear()
            self.setCurrentFileName()

    def fileOpen(self):
            fname = QFileDialog.getOpenFileName(self, "Open File...", exedir,
                    "Input Files (*.inp);; All Files (*)",QFileDialog.DontUseNativeDialog)
                    #None,filter,QFileDialog.DontUseNativeDialog)
            if fname:
                self.load(fname)

    def load(self,f):

        if not QFile.exists(f):
            return False

        fh = QFile(f)
        if not fh.open(QFile.ReadOnly):
            return False

        instr = QTextStream(fh)
        # FREEZE CURSOR
        QApplication.setOverrideCursor(Qt.WaitCursor)

        self.inputTextEdit.setPlainText(instr.readAll())

        # NORMAL CURSOR
        QApplication.restoreOverrideCursor()

        self.setCurrentFileName(f)
        #self.tabs.setCurrentIndex(2)
        self.action_SaveAs.setEnabled(True)
        self.chargeSpinBox.setValue(0)
        self.statusBar().showMessage("File loaded", 2000)
        return True

    def maybeSave(self):
        if not self.inputTextEdit.document().isModified():
            return True

        if self.fileName.startswith(':/'):
            return True

        ret = QMessageBox.warning(self, "Application",
                "The input has been modified.\n"
                "Do you want to save your changes?",
                QMessageBox.Save | QMessageBox.Discard |
                        QMessageBox.Cancel)

        if ret == QMessageBox.Save:
            return self.fileSave()

        if ret == QMessageBox.Cancel:
            return False

        return True

    def maybeRun(self):
        if self.fileName.startswith(':/'):
            self.tabs.setCurrentIndex(2)
            if self.inputTextEdit.document().isModified():
                ret = QMessageBox.warning(self, "Application",
                    "The input has been modified.\n"
                    "Do you want to save your changes?",
                    QMessageBox.Save | QMessageBox.Discard |
                            QMessageBox.Cancel)

            else:
                self.noinputerror()
            return False

        if self.inputTextEdit.document().isModified():
            ret = QMessageBox.warning(self, "Application",
                    "The input has been modified.\n"
                    "Do you want to save your changes?",
                    QMessageBox.Save | QMessageBox.Discard |
                            QMessageBox.Cancel)

            if ret == QMessageBox.Save:
                self.fileSave()
                return True

            if ret == QMessageBox.Discard:
                return True

            if ret == QMessageBox.Cancel:
                return False
        else:
            return True

    def testinp(self):
        self.tabs.setCurrentIndex(2)

        if not self.fileName.startswith(':/'):
            from get_basis import Z_tablelist
            Z_table = Z_tablelist()
            f = open(self.fileName, 'r')
            lines = f.readlines()
            errorinp = 0
            numoflines = 0
            for line in lines:
                lines = line.strip()
                if lines:
                    try:
                        if len(line.split()) == 4 and line.split()[0].upper() in Z_table and line.split()[1:4] and all(isinstance(item, (float,int)) for item in map(float,line.split()[1:4])):
                            numoflines += 1

                        else:
                            print(line.split(),len(line.split()),line.split()[0].upper(),line.split()[1:4])
                            errorinp += 1

                        #print(line.split())
                        #print (len(line.split()), "columns")
                    except ValueError:
                        errorinp += 1


            f.close()
            if errorinp == 0 and numoflines > 0:
                return True
            else:
                return False

    def fileSave(self):

        if self.fileName.startswith(':/'):
            return self.fileSaveAs()

        fh = QFile(self.fileName)

        if not self.fileName:
            return self.fileSaveAs()

        if not fh.open(QFile.WriteOnly | QFile.Text):
            QMessageBox.warning(self, "Recent Files",
            "Cannot write file %s:\n%s." % (self.fileName, fh.errorString()))
            return

        outstr = QTextStream(fh)

        QApplication.setOverrideCursor(Qt.WaitCursor)

        outstr << self.inputTextEdit.toPlainText()

        QApplication.restoreOverrideCursor()

        self.inputTextEdit.document().setModified(False)

        self.statusBar().showMessage("File saved", 2000)

    def fileSaveAs(self):
        fn = QFileDialog.getSaveFileName(self, "Save as...", exedir+"/"+self.fileName.split('/')[-1].split('.')[0],
                "Input File (*.inp);;All Files (*)")

        if not fn:
            return False

        lfn = fn.lower()
        if not lfn.endswith(('.inp')):
            # The default.
            fn += '.inp'

        self.setCurrentFileName(fn)
        return self.fileSave()


    def fileSaveOutputAs(self):
        #self.outputTextEdit.setReadOnly(False)
        fn = QFileDialog.getSaveFileName(self, "Save as...", exedir+"/"+self.fileName.split('/')[-1].split('.')[0],
                "Output File (*.out);;All Files (*)")

        if not fn:
            return False

        lfn = fn.lower()
        if not lfn.endswith(('.out')):
            # The default.
            fn += '.out'
        self.f = open(fn,'w')
        return self.fileSaveOutput()


    def fileSaveOutput(self):

        filedata = self.outputTextEdit.toPlainText()
        self.f.write(filedata)
        self.f.close()

        self.statusBar().showMessage("File saved", 2000)
        
    def MODensityCall(self):
        
        # FREEZE CURSOR
        QApplication.setOverrideCursor(Qt.WaitCursor)

        if self.MOquestion():
            self.fnCube = QFileDialog.getSaveFileName(self, "Save Cube File as...",exedir+"/"+self.fileName.split('/')[-1].split('.')[0],
                                                      "Cube File (*.cube);;All Files (*)",QFileDialog.DontUseNativeDialog)

            if not self.fnCube:
                return False
            else:
                lfn = self.fnCube.lower()
                if not lfn.endswith(('.cube')):
                    # The default.
                    self.fnCube += '.cube'
                    
                if sys.platform.startswith('win32'):
                    self.libCubePath = exedir+"\lib\_CUBE.dll"
                elif sys.platform.startswith('linux'):
                    self.libCubePath = exedir+"/lib/_CUBE.so"
                elif sys.platform.startswith('darwin'):
                    self.libCubePath = exedir+"/lib/_CUBE.dylib"

                
        
                from cubeCall import cubeCall
                cubeCall(self.R,C,P,self.A,self.D,self.Center,self.CartAng,self.lenL,self.locL,self.Z,self.fnCube,self.N,self.n,self.MOstart,self.MOend,self.plot_type,self.libCubePath,logger)
        # NORMAL CURSOR
        QApplication.restoreOverrideCursor()

    def openRecentFile(self):
        action = self.sender()
        if action:
            self.load(action.data())
            
            
    def updateRecentFileActions(self):
        settings = QSettings(PACKAGENAME, 'Recent Files')
        files = settings.value('recentFileList', [])

        numRecentFiles = min(len(files), self.MaxRecentFiles)

        for i in range(numRecentFiles):
            text = "&%d %s" % (i + 1, self.strippedName(files[i]))
            self.recentFileActsInp[i].setText(text)
            self.recentFileActsInp[i].setData(files[i])
            self.recentFileActsInp[i].setVisible(True)

        for j in range(numRecentFiles, self.MaxRecentFiles):
            self.recentFileActsInp[j].setVisible(False)

        self.separatorAct.setVisible((numRecentFiles > 0))

    def strippedName(self, fullFileName):
        return QFileInfo(fullFileName).fileName()

    def filePrint(self):
        printer = QPrinter(QPrinter.HighResolution)
        dlg = QPrintDialog(printer, self)

        if self.outputTextEdit.textCursor().hasSelection():
            dlg.addEnabledOption(QAbstractPrintDialog.PrintSelection)

        dlg.setWindowTitle("Print Document")

        if dlg.exec_() == QDialog.Accepted:
            self.outputTextEdit.print_(printer)

        del dlg

    def filePrintPreview(self):
        printer = QPrinter(QPrinter.HighResolution)
        preview = QPrintPreviewDialog(printer, self)
        preview.paintRequested.connect(self.printPreview)
        preview.exec_()

    def printPreview(self, printer):
        self.outputTextEdit.print_(printer)

    def about(self, checked=None):
        dialog = design.aboutDialog(self)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.exec_()

    #class aboutwin(Q)


    def about2(self):
        #self.ui = aboutDialog()
        if self.aboutDialog is None:
            self.aboutDialog = design.aboutDialog(self)

        if self.aboutDialog.exec_():
            settings = QSettings(self.design.aboutDialog.format())


        #if self.locationDialog is None:
        #    self.locationDialog = LocationDialog(self)

        #if self.locationDialog.exec_():
        #    settings = QSettings(self.locationDialog.format())

    def about2(self):
        QMessageBox.about(self, "About","""
                <HTML>
                <p><font size = "5"><center><b>Simple Hartree-Fock Program Calculation To Analyse Levels Energy
                                And Orbital Cube View of Atoms and Small Molecules.</b></center></p>
                                <BR>
                                <BR>
                                <BR>
                                <p><center> Weverson R. Gomes (gomeswr@gmail.com)</center></p>
                                """)
    def closeEvent(self, event):

        if self.maybeSave():
            event.accept()
            close('all')
        else:
            event.ignore()

    def inputerror(self):
        QMessageBox.critical(self, "ERROR",
        "Input ERROR.",QMessageBox.Ok)
        return

    def noinputerror(self):
        QMessageBox.critical(self, "NO INPUT",
        "You didn't enter any input.",QMessageBox.Ok)
        return

    def plotorbital(self):
        s = {1:[5],2:[4,6],3:[4,5,6],4:[3,4,6,7],5:[3,4,5,6,7]}
        err = 0.01

        eround = [ round(elem, 5) for elem in e ]
        dicte = Counter(eround)

        x = []
        y = [repvalue for evalues,deg in sorted(dicte.items()) for repvalue in repeat(evalues,deg)]

        for evalues,deg in sorted(dicte.items()):
            x.extend(s[deg])

        rcParams['interactive']  = True
        f1 = figure(figsize=(5,4))
        ax1 = f1.add_subplot(111)
        ax1.axes.get_xaxis().set_visible(False)
        ylabel('Electronic Energy (a.u)',fontsize = 18)
        xlim(-2,12)
        f1.suptitle(self.fileName.split('/')[-1])
        ax1.plot(x,y,'_',markersize=20,color='k',mew=1.2,label=self.basis+'/'+'charge='+str(self.charge))
#plt.plot(self.countplot,self.ENplot,'o-',label="label")
        legend()

    def plotconvergence(self):
        #ENplot =  ENlist[:]
        #countplot = countlist[:]
        rcParams['interactive']  = True
        f2 = figure(figsize=(5,4))
        f2.suptitle(self.fileName.split('/')[-1])
        ax2 = f2.add_subplot(111)
        #ax2.hold(True)
        ax2.plot(countlist,ENlist,'*-',label=self.basis+'/'+'charge='+str(self.charge))
        xlabel('Steps',fontsize = 18)
        ylabel('Electronic Energy (a.u)',fontsize = 18)
        legend()


class MovieSplashScreen(QSplashScreen):

    def __init__(self, movie, parent = None):

        movie.jumpToFrame(0)
        pixmap = QPixmap(movie.frameRect().size())

        QSplashScreen.__init__(self, pixmap)
        self.movie = movie
        self.movie.frameChanged.connect(self.repaint)

    def showEvent(self, event):
        self.movie.start()

    def hideEvent(self, event):
        self.movie.stop()

    def paintEvent(self, event):

        painter = QPainter(self)
        pixmap = self.movie.currentPixmap()
        self.setMask(pixmap.mask())
        painter.drawPixmap(0, 0, pixmap)

    def sizeHint(self):

        return self.movie.scaledSize()


class aboutDialog2(QDialog):
    def __init__(self, parent=None):
        super(aboutDialog, self).__init__(parent)

        aboutText = """
                <HTML>
                <p><font size = "5"><center><b>Easy Hartree-Fock Program Calculation</b></center></p>
                <p><font size = "5"><b><center>To Levels Energy Analyse And Orbital Cube View.</b></center></p>
                                 <BR>
                                 <BR>
                                 <BR>
               <p><center> Weverson R. Gomes (gomeswr@gmail.com)</center></p>
                                 """


        self.lpix = QLabel(self)
        self.label = QLabel(self)
        self.label.setGeometry(QRect(30, 100, 341, 61))
        self.label.setObjectName("label")
        self.label.setText(aboutText)

        self.licenseButton = QPushButton(self)
        #self.licenseButton.setGeometry(QRect(10, 150, 301, 41))
        self.licenseButton.setObjectName("licenseButton")
        self.licenseButton.setText("&License")

        self.closeButton = QPushButton(self)
        #self.pushButton.setGeometry(QRect(10, 150, 301, 41))
        self.closeButton.setObjectName("pushButton")
        self.closeButton.setText("&Close")
        #QObject.connect(self.closeButton, SIGNAL("clicked()"), self.closeabout)
        #self.closeButton.clicked.connect(self.closeabout())
        self.connect(self.closeButton, SIGNAL("clicked()"), self.aboutClose)

        #hlayout = QFormLayout()
        #hlayout.addRow(self.pushButton)
        #hlayout.addRow(self.licenseButton)
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.licenseButton)
        hlayout.addWidget(self.closeButton)


        layout = QGridLayout()
        layout.addWidget(self.lpix,0,0,1,3, Qt.AlignCenter)
        layout.addWidget(self.label,1,0,1,3)
        layout.addLayout(hlayout,2,2)
        self.setLayout(layout)

        self.pixmap = QPixmap("fora.png")
        self.lpix.setPixmap(self.pixmap)

        QMetaObject.connectSlotsByName(self)


def main():
    app = QApplication(sys.argv)

    #splash_pix = QPixmap('./images/fora.png')
    #splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    print(exedir+"/images/splash/ezgif.com-gif-maker_2.gif")
    movie = QMovie(exedir+"/images/splash/ezgif.com-gif-maker_2.gif")
    splash = MovieSplashScreen(movie)
    #splash.setMask(splash_pix.mask())
    splash.show()

    starttime = time.time()
    while time.time() < starttime + 4:
        app.processEvents()


    form = sashaprogram()
    form.show()
    splash.finish(form)
    #app.exec_()
    sys.exit(app.exec_())


if __name__ == '__main__':


    for i in range(i0, 1):
        progressBar.setValue(i)
        t = time.time()
        while time.time() < t + 0.05:
            app.processEvents()

    # Simulate something that takes time
    form = sashaprogram()
    form.show()
    splash.finish(form)
    #app.exec_()
    sys.exit(app.exec_())
