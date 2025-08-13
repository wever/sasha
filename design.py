# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'shasha_esse2.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!


import sip
sip.setapi('QString', 2)

from PyQt4.QtGui import  QApplication,QDesktopWidget,QWidget,QFont,QDialog,QDialogButtonBox,QIcon,QVBoxLayout,QHBoxLayout,QTabWidget,QGridLayout,QLabel,QSizePolicy,QSpacerItem,QComboBox,QSpinBox,QCheckBox,QTextEdit,QPushButton,QMenuBar,QMenu,QStatusBar,QAction,QKeySequence,QPixmap,QTextBrowser,qApp
from PyQt4.QtCore import QRect,Qt,QMetaObject,SIGNAL,SLOT,QObject
import sys
from get_basis import basislist
import sasha_rc
import os

def _fromUtf8(s):
    return s

try:
    _encoding = QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig)


# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    exedir= os.path.dirname(sys.executable)
elif __file__:
    exedir = os.path.dirname(__file__)

if sys.platform.startswith('darwin'):
    rsrcPath = exedir+"/images/mac"
elif sys.platform.startswith('linux'):
    rsrcPath = exedir+"/images/linux"
else:
    rsrcPath = exedir+"/images/win"


class Ui_MainWindow(object):

    def center(self):

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def setupUi(self,MainWindow):
        self.setWindowIcon(QIcon(exedir+'/images/logo.png'))
        MainWindow.setObjectName(_fromUtf8("SASHA PACKAGE"))
        MainWindow.resize(800, 600)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout = QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.tabs = QTabWidget(self.centralwidget)
        self.tabs.setObjectName(_fromUtf8("tabs"))
        self.tab = QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.gridLayout = QGridLayout(self.tab)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.methodLabel = QLabel(self.tab)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.methodLabel.sizePolicy().hasHeightForWidth())
        self.methodLabel.setSizePolicy(sizePolicy)
        font = QFont()
        font.setPointSize(10)
        self.methodLabel.setFont(font)
        self.methodLabel.setObjectName(_fromUtf8("methodLabel"))
        self.gridLayout.addWidget(self.methodLabel, 1, 0, 1, 1)
        spacerItem = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 6, 0, 1, 1)
        spacerItem1 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 0, 1, 1, 1)
        self.comboBox = QComboBox(self.tab)
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.comboBox.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.comboBox, 2, 0, 1, 1)
        self.basissetLabel = QLabel(self.tab)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.basissetLabel.sizePolicy().hasHeightForWidth())
        self.basissetLabel.setSizePolicy(sizePolicy)
        font = QFont()
        font.setPointSize(10)
        self.basissetLabel.setFont(font)
        self.basissetLabel.setObjectName(_fromUtf8("basissetLabel"))
        self.gridLayout.addWidget(self.basissetLabel, 1, 1, 1, 1)
        spacerItem2 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem2, 3, 0, 1, 1)
        spacerItem3 = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem3, 0, 0, 1, 1)
        self.comboBox_2 = QComboBox(self.tab)
        self.comboBox_2.setObjectName(_fromUtf8("comboBox_2"))

        global basis
        basis = basislist()

        for i in range(len(basis)):
            self.comboBox_2.addItem(_fromUtf8(""))

        self.gridLayout.addWidget(self.comboBox_2, 2, 1, 1, 1)
        self.unitsComboBox = QComboBox(self.tab)
        self.unitsComboBox.setObjectName(_fromUtf8("unitsComboBox"))
        self.gridLayout.addWidget(self.unitsComboBox, 5, 0, 1, 1)
        self.unitsComboBox.addItem("Angstroms")
        self.unitsComboBox.addItem("Bohr")
        self.unitsComboBox.setCurrentIndex(0)
        self.unitsLabel = QLabel(self.tab)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.unitsLabel.sizePolicy().hasHeightForWidth())
        self.unitsLabel.setSizePolicy(sizePolicy)
        font = QFont()
        font.setPointSize(10)
        self.unitsLabel.setFont(font)
        self.unitsLabel.setObjectName(_fromUtf8("unitsLabel"))
        self.gridLayout.addWidget(self.unitsLabel, 4, 0, 1, 1)

        self.chargeLabel = QLabel(self.tab)
        self.chargeSpinBox = QSpinBox(self.tab)
        self.chargeSpinBox.setObjectName(_fromUtf8("chargeSpinBox"))
        self.gridLayout.addWidget(self.chargeSpinBox, 5, 1, 1, 1)
        font = QFont()
        font.setPointSize(10)
        self.chargeLabel.setFont(font)
        self.chargeLabel.setObjectName(_fromUtf8("unitsLabel"))
        self.gridLayout.addWidget(self.chargeLabel, 4, 1, 1, 1)
        self.chargeSpinBox.setMinimum(-10)
        self.chargeSpinBox.setMaximum(10)
        self.chargeSpinBox.setValue(0)

        self.tabs.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.gridLayout_2 = QGridLayout(self.tab_2)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.checkBox = QCheckBox(self.tab_2)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.gridLayout_2.addWidget(self.checkBox, 0, 0, 1, 1)
        #self.gridLayout_2.addWidget(self.checkBox, 1, 0, 1, 1)
        self.checkBox_2 = QCheckBox(self.tab_2)
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.gridLayout_2.addWidget(self.checkBox_2, 1, 0, 1, 1)
        #self.gridLayout_2.addWidget(self.checkBox_2, 1, 2, 1, 1)
        self.checkBox_3 = QCheckBox(self.tab_2)
        self.checkBox_3.setObjectName(_fromUtf8("checkBox_3"))
        self.gridLayout_2.addWidget(self.checkBox_3, 2, 0, 1, 1)
        #self.gridLayout_2.addWidget(self.checkBox_3, 0, 0, 1, 1)
        self.checkBox_4 = QCheckBox(self.tab_2)
        self.checkBox_4.setObjectName(_fromUtf8("checkBox_4"))
        self.gridLayout_2.addWidget(self.checkBox_4, 0, 2, 1, 1)
        self.checkBox_5 = QCheckBox(self.tab_2)
        self.checkBox_5.setObjectName(_fromUtf8("checkBox_5"))
        self.gridLayout_2.addWidget(self.checkBox_5, 1, 2, 1, 1)
        #self.gridLayout_2.addWidget(self.checkBox_5, 2, 0, 1, 1)
        self.checkBox_6 = QCheckBox(self.tab_2)
        self.checkBox_6.setObjectName(_fromUtf8("checkBox_6"))
        self.gridLayout_2.addWidget(self.checkBox_6, 2, 2, 1, 1)
        self.tabs.addTab(self.tab_2, _fromUtf8(""))

        self.tab_5 = QWidget()
        self.tab_5.setObjectName(_fromUtf8("tab_5"))
        self.verticalLayout_3 = QVBoxLayout(self.tab_5)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        #self.inputTextEdit = QPlainTextEdit(self.tab_5)
        self.inputTextEdit = QTextEdit(self.tab_5)
        #self.inputTextEdit = QTextEdit(self.tab_5)
        #self.inputTextEdit = QTextBrowser(self.tab_5)
        self.inputTextEdit.setObjectName(_fromUtf8("inputTextEdit"))
        self.inputTextEdit.setFontFamily("Courier")
        #self.inputTextEdit.setFontFamily("monospace")
        self.inputTextEdit.setFontPointSize(12)
        self.verticalLayout_3.addWidget(self.inputTextEdit)
        self.tabs.addTab(self.tab_5, _fromUtf8(""))

        self.tab_3 = QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.verticalLayout_2 = QVBoxLayout(self.tab_3)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        #self.outputTextEdit = QPlainTextEdit(self.tab_3)
        self.outputTextEdit = QTextEdit(self.tab_3)
        #self.outputTextEdit.setAlignment(Qt.AlignCenter)
        #self.outputTextEdit.setCurrentFont(QFont('fixed'))
        #self.outputTextEdit.setStyleSheet('font-size: 10pt; font-family: Consolas;')
        #self.outputTextEdit = QTextBrowser(self.tab_3)
        self.outputTextEdit.setFontFamily("Courier")
        #self.outputTextEdit.setFontFamily("monospace")
        self.outputTextEdit.setFontPointSize(12)
        self.outputTextEdit.setReadOnly(True)
        self.outputTextEdit.setObjectName(_fromUtf8("outputTextEdit"))
        self.verticalLayout_2.addWidget(self.outputTextEdit)
        self.tabs.addTab(self.tab_3, _fromUtf8(""))
        self.tab_4 = QWidget()
        self.tab_4.setObjectName(_fromUtf8("tab_4"))
        self.gridLayout_3 = QGridLayout(self.tab_4)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.orbitalButton = QPushButton(self.tab_4)
        self.orbitalButton.setObjectName(_fromUtf8("orbitalButton"))
        self.gridLayout_3.addWidget(self.orbitalButton, 0, 0, 1, 1)
        self.convergenceButton = QPushButton(self.tab_4)
        self.convergenceButton.setObjectName(_fromUtf8("convergenceButton"))
        self.gridLayout_3.addWidget(self.convergenceButton, 0, 1, 1, 1)
        self.MODensityButton = QPushButton(self.tab_4)
        self.MODensityButton.setObjectName(_fromUtf8("MODensityButton"))
        self.gridLayout_3.addWidget(self.MODensityButton, 1, 0, 1, 1)
        self.printoutputButton = QPushButton(self.tab_4)
        self.printoutputButton.setObjectName(_fromUtf8("printoutputButton"))
        self.gridLayout_3.addWidget(self.printoutputButton, 1, 1, 1, 1)

        self.tabs.addTab(self.tab_4, _fromUtf8(""))
        self.verticalLayout.addWidget(self.tabs)

        self.orbitalButton.setEnabled(False)
        self.convergenceButton.setEnabled(False)
        self.MODensityButton.setEnabled(False)
        self.printoutputButton.setEnabled(False)
        
        self.center()


        # RUN BUTTON AND ACTION RUN SASHA PACKAGE
        self.runButton = QPushButton(self.centralwidget)#,shortcut=Qt.CTRL + Qt.Key_B)
        self.runButton.setObjectName(_fromUtf8("runButton"))
        self.runButton.setEnabled(True)
        self.killButton = QPushButton(self.centralwidget)
        self.killButton.setObjectName(_fromUtf8("killButton"))
        self.killButton.setEnabled(False)
    
        self.verticalLayout.addWidget(self.runButton)
        self.verticalLayout.addWidget(self.killButton)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setGeometry(QRect(0, 0, 800, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menu_FIle = QMenu(self.menubar)
        self.menu_FIle.setObjectName(_fromUtf8("menu_FIle"))
        self.menu_Help = QMenu(self.menubar)
        self.menu_Help.setObjectName(_fromUtf8("menu_Help"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)


        # OPEN BUTTON AND ACTION OPEN WINDOW
        self.action_Open = QAction(
                QIcon.fromTheme('document-open',
                        QIcon(rsrcPath + '/fileopen.png')),
                "&Open", self, shortcut=QKeySequence.Open)
        self.menu_FIle.addAction(self.action_Open)
        #self.action_Open = QAction(QIcon('open.png'), 'Open', self)
        #self.action_Open.setShortcut('Ctrl+O')
        #self.action_Open.setStatusTip('Open new File')
        #self.action_Open.setObjectName(_fromUtf8("action_Open"))

       # SAVE BUTTON AND ACTION SAVE INPUT
        self.action_Save = QAction(
                QIcon.fromTheme('document-save',
                        QIcon(rsrcPath + '/filesave.png')),
                "&Save Input", self, shortcut=QKeySequence.Save,
                enabled=False)
        self.menu_FIle.addAction(self.action_Save)
  
       # SAVE BUTTON AND ACTION SAVE INPUT AS
        self.action_SaveAs = QAction("Save &Input As...", self,
                priority=QAction.LowPriority,
                shortcut=Qt.CTRL + Qt.SHIFT + Qt.Key_S,
                enabled=False)
        self.menu_FIle.addAction(self.action_SaveAs)
        self.menu_FIle.addSeparator()

       # SAVE BUTTON AND ACTION SAVE OUTPUT AS
        self.action_SaveOutputAs = QAction("Save Output &As...", self,
                priority=QAction.LowPriority,
                shortcut=Qt.CTRL + Qt.Key_W,
                enabled=False)
        self.menu_FIle.addAction(self.action_SaveOutputAs)
        
       # RECENT FILES ACTION 
        self.separatorAct = self.menu_FIle.addSeparator()

        self.MaxRecentFiles = 5
        self.recentFileActsInp = []

        for i in range(self.MaxRecentFiles):
            self.recentFileActsInp.append(
                    QAction(self, visible=False,
                            triggered=self.openRecentFile))
 
 
        for i in range(self.MaxRecentFiles):
            self.menu_FIle.addAction(self.recentFileActsInp[i])
     
        self.menu_FIle.addSeparator()

        # PRINT AND PREVIEW BUTTON
        self.action_Print = QAction(
                QIcon.fromTheme('document-print',
                        QIcon(rsrcPath + '/fileprint.png')),
                "&Print...", self, priority=QAction.LowPriority,
                shortcut=QKeySequence.Print)
        self.menu_FIle.addAction(self.action_Print)
        self.menu_FIle.addAction(self.action_Print)
  
        self.action_PrintPreview = QAction(
                QIcon.fromTheme('fileprint',
                        QIcon(rsrcPath + '/fileprint.png')),
                "Print Preview...", self,
                shortcut=Qt.CTRL + Qt.SHIFT + Qt.Key_P)
        self.menu_FIle.addAction(self.action_PrintPreview)
        self.menu_FIle.addSeparator()

        # EXIT BUTTON AND ACTION ASK FOR QUIT
        self.actionExit = QAction("&Exit", self,
                shortcut=QKeySequence.Quit)
        self.menu_FIle.addAction(self.actionExit)

        # BUTTON ABOUT
        self.action_About = QAction(MainWindow)
        self.action_AboutQt = QAction(MainWindow)

        self.action_About.setObjectName(_fromUtf8("action_About"))
        self.action_AboutQt.setObjectName(_fromUtf8("action_AboutQt"))
        #self.menu_FIle.addAction(self.action_Open)
        #self.menu_FIle.addSeparator()
        #self.menu_FIle.addAction(self.action_Save)
        #self.menu_FIle.addAction(self.action_SaveAs)
        #self.menu_FIle.addSeparator()
        #self.menu_FIle.addAction(self.actionExit)
        self.menu_Help.addAction(self.action_About)
        self.menu_Help.addAction(self.action_AboutQt)
        self.menubar.addAction(self.menu_FIle.menuAction())
        self.menubar.addAction(self.menu_Help.menuAction())

        self.retranslateUi(MainWindow)
        QObject.connect(self.actionExit, SIGNAL(_fromUtf8("activated()")), MainWindow.close)
        QMetaObject.connectSlotsByName(MainWindow)
        


    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("SaSHA PACKAGE", "SaSHA PACKAGE", None))
        self.methodLabel.setText(_translate("MainWindow", "Method:", None))
        self.chargeLabel.setText(_translate("MainWindow", "Charge:", None))
        self.comboBox.setItemText(0, _translate("MainWindow", "Hartree-Fock", None))
        self.basissetLabel.setText(_translate("MainWindow", "Basis Set:", None))


        for i,bs in enumerate(basis):
            self.comboBox_2.setItemText(i, _translate("MainWindow", bs, None))

        self.comboBox_2.setCurrentIndex(0)
    
        self.unitsLabel.setText(_translate("MainWindow", "Input Units:", None))
        self.tabs.setTabText(self.tabs.indexOf(self.tab), _translate("MainWindow", "MENU", None))
        self.checkBox.setText(_translate("MainWindow", "Print Overlap Matrix", None))
        self.checkBox_2.setText(_translate("MainWindow", "Print Kinetic Energy Matrix", None))
        self.checkBox_3.setText(_translate("MainWindow", "Print Three-Center Nuclear Attraction Matrix", None))
        self.checkBox_4.setText(_translate("MainWindow", "Print Two-Electron Repulsion Matrix", None))
        self.checkBox_5.setText(_translate("MainWindow", "Save Output", None))
        self.checkBox_6.setText(_translate("MainWindow", "Print Molecular Orbital / Electron Density (Cube Format)", None))

        # TABS
        self.tabs.setTabText(self.tabs.indexOf(self.tab_2), _translate("MainWindow", "PRINT TO FILE", None))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_3), _translate("MainWindow", "OUTPUT", None))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_4), _translate("MainWindow", "ANALYSIS", None))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_5), _translate("MainWindow", "INPUT", None))
        self.tabs.setCurrentIndex(0)

        self.menu_FIle.setTitle(_translate("MainWindow", "&File", None))
        self.menu_Help.setTitle(_translate("MainWindow", "&Help", None))

        # ACTION BUTTON GRAPH ORBITALS
        self.convergenceButton.setText(_translate("MainWindow", "Convergence", None))
        self.convergenceButton.setStatusTip(self.tr("Show Convergence Steps"))

        # ACTION BUTTON GRAPH ORBITALS
        self.orbitalButton.setText(_translate("MainWindow", "Orbitals Energies", None))
        self.orbitalButton.setStatusTip(self.tr("Show Orbital Energies"))

        # ACTION BUTTON GRAPH MO / DENSITY
        self.MODensityButton.setText(_translate("MainWindow", "Molecular Orbital / Electron Density", None))
        self.MODensityButton.setStatusTip(self.tr("Print Molecular Orbital / Electron Density"))
        
        # ACTION BUTTON SAVE OUTPUT
        self.printoutputButton.setText(_translate("MainWindow", "Save Output", None))
        self.printoutputButton.setStatusTip(self.tr("Save Output log"))

        # ACTION BUTTON RUN
        self.runButton.setText(_translate("MainWindow", "&RUN", None))

        # ACTION BUTTON KILL
        self.killButton.setText(_translate("MainWindow", "&KILL", None))

        # ACTION BUTTON OPEN
        self.action_Open.setText(_translate("MainWindow", "&Open", None))
        self.action_Open.setStatusTip(self.tr("Open Input File"))

        # ACTION BUTTON SAVE
        self.action_Save.setText(_translate("MainWindow", "&Save Input", None))
        self.action_Save.setStatusTip(self.tr("Save Input File"))

        # ACTION BUTTON SAVE AS
        self.action_SaveAs.setText(_translate("MainWindow", "Save &Input As", None))
        self.action_SaveAs.setStatusTip(self.tr("Save Input File As"))

        # ACTION BUTTON SAVE OUTPUT AS
        self.action_SaveOutputAs.setText(_translate("MainWindow", "Save Output &As", None))
        self.action_SaveOutputAs.setStatusTip(self.tr("Save Output File As"))

        # ACTION PRINT AND PREVIEW
        self.action_Print.setText(_translate("MainWindow", "&Print Output", None))
        self.action_Print.setStatusTip(self.tr("Print Output File"))
        self.action_PrintPreview.setText(_translate("MainWindow", "Preview &Output", None))
        self.action_PrintPreview.setStatusTip(self.tr("Preview Output File"))

        # ACTION BUTTON EXIT
        self.actionExit.setText(_translate("MainWindow", "&Exit", None))
        self.actionExit.setStatusTip(self.tr("Quit Program"))

        # ACTION BUTTON ABOUT
        self.action_About.setText(_translate("MainWindow", "&About", None))
        self.action_About.setStatusTip(self.tr("About Program"))

        # ACTION BUTTON ABOUT QT
        self.action_AboutQt.setText(_translate("MainWindow", "About &Qt", None))
        self.action_AboutQt.setStatusTip(self.tr("Show the Qt library's About box"))
        self.connect(self.action_AboutQt, SIGNAL("triggered()"), qApp, SLOT("aboutQt()"))

class MODialog(QDialog):
    def __init__(self, N, n_electrons, parent=None):
        super(MODialog, self).__init__(parent)
        self.N = N
        self.setObjectName(_fromUtf8("M.O and E.D CUBE"))
        self.resize(400, 123)
        self.gridLayout_MO2 = QGridLayout(self)
        self.gridLayout_MO2.setObjectName(_fromUtf8("gridLayout_MO2"))
        self.comboBox_MO_DENSITY = QComboBox(self)
        self.comboBox_MO_DENSITY.setObjectName(_fromUtf8("comboBox_MO_DENSITY"))
        self.gridLayout_MO2.addWidget(self.comboBox_MO_DENSITY, 1, 1, 1, 1)
        self.comboBox_end = QComboBox(self)
        self.comboBox_end.setObjectName(_fromUtf8("comboBox_end"))
        self.gridLayout_MO2.addWidget(self.comboBox_end, 5, 2, 1, 1)
        spacerItemMO = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout_MO2.addItem(spacerItemMO, 6, 0, 1, 1)
        spacerItem1MO = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout_MO2.addItem(spacerItem1MO, 0, 1, 1, 1)
        self.buttonBoxOK = QDialogButtonBox(self)
        self.buttonBoxOK.setOrientation(Qt.Horizontal)
        self.ret = self.buttonBoxOK.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBoxOK.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout_MO2.addWidget(self.buttonBoxOK, 9, 0, 1, 3)
        self.comboBox_start = QComboBox(self)
        self.comboBox_start.setObjectName(_fromUtf8("comboBox_start"))
        self.gridLayout_MO2.addWidget(self.comboBox_start, 5, 0, 1, 1)
        spacerItem2MO = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.gridLayout_MO2.addItem(spacerItem2MO, 6, 2, 1, 1)
        self.label_MO3 = QLabel(self)
        self.label_MO3.setObjectName(_fromUtf8("label_MO3"))
        self.gridLayout_MO2.addWidget(self.label_MO3, 4, 2, 1, 1)
        self.label_MO2 = QLabel(self)
        self.label_MO2.setObjectName(_fromUtf8("label_MO2"))
        self.gridLayout_MO2.addWidget(self.label_MO2, 4, 0, 1, 1)
        self.setFixedSize(400, 150)
        self.comboBox_MO_DENSITY.addItem("Molecular Orbital Cube")
        self.comboBox_MO_DENSITY.addItem("Electron Density Cube")
        self.label_MO2.setText("From Orbital:")
        self.label_MO3.setText("To Orbital:")

        print(n_electrons)
        orb_occ = n_electrons//2
        orb_virtual = N - orb_occ
        for i in range(1,self.N+1):
            if i <= orb_occ:
                if i == orb_occ:
                    printformatted = "%s %s %s %s" % ("HOMO"," (",str(i),")")
                    self.comboBox_start.addItem(printformatted)
                    self.comboBox_end.addItem(printformatted)
                else:
                    printformatted = "%s %+d %s %s %s" % ("HOMO",i-orb_occ," (",str(i),")")
                    #self.comboBox_start.addItem("HOMO"+ str(i-orb_occ)+ " ("+str(i)+")")
                    self.comboBox_start.addItem(printformatted)
                    self.comboBox_end.addItem("HOMO"+ str(i-orb_occ)+ " ("+str(i)+")")
            else:
                if i == orb_occ+1:
                    self.comboBox_start.addItem("LUMO"+" ("+str(i)+")")
                    self.comboBox_end.addItem("LUMO"+" ("+str(i)+")")
                else:
                    self.comboBox_start.addItem("LUMO"+ str(i-orb_occ)+ " ("+str(i)+")")
                    self.comboBox_end.addItem("LUMO"+ str(i-orb_occ)+ " ("+str(i)+")")

        self.connect(self.buttonBoxOK, SIGNAL(_fromUtf8("accepted()")), self.accept)
        self.connect(self.buttonBoxOK, SIGNAL(_fromUtf8("rejected()")), self.reject)
        self.comboBox_start.currentIndexChanged.connect(self.changecombo)
        self.comboBox_MO_DENSITY.currentIndexChanged.connect(self.changecomboMO)

        
    def changecombo(self):
        self.comboBox_end.clear()
        for j in range(self.comboBox_start.currentIndex(),self.N):
            self.comboBox_end.addItem(str(j))

    def changecomboMO(self):
        if self.comboBox_MO_DENSITY.currentText() == "Electron Density Cube":
            self.comboBox_start.setEnabled(False)
            self.comboBox_end.setEnabled(False)
        else:
            self.comboBox_start.setEnabled(True)
            self.comboBox_end.setEnabled(True)


class aboutDialog(QDialog):
    def __init__(self,parent=None):
        super(aboutDialog, self).__init__(parent)
        
        aboutText = """ 
                <HTML>
                <p><font size = "4"><center><b>Simple Hartree-Fock Program Calculation</b></center></p> 
                <p><font size = "4"><b><center>To Energy Levels Analyse And</b></center></p>
                <p><font size = "4"><b><center>Orbital Cube View of Atoms and Small Molecules.</b></center></p>
                                 <BR>
                                 <BR>
                                 <BR>
               <p><center> Weverson R. Gomes (gomeswr@gmail.com)</center></p>
                                 """
        self.setWindowTitle("About Sasha")
        self.resize(400, 301)


        self.lpix = QLabel(self)

        self.label = QLabel(self)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.label.setText(aboutText)

        self.labelblankline = QLabel(self)
        self.labelblankline.setGeometry(QRect(30, 100, 341, 61))
        self.labelblankline.setObjectName("labelblankline")


        self.licenseButton = QPushButton(self)
        #self.licenseButton.setGeometry(QRect(10, 150, 301, 41))
        self.licenseButton.setObjectName("licenseButton")
        self.licenseButton.setText("&License")
        
        self.closeButton = QPushButton(self)
        #self.pushButton.setGeometry(QRect(10, 150, 301, 41))
        self.closeButton.setObjectName("pushButton")
        self.closeButton.setText("&Close")
        #self.closeButton.setIcon(QIcon(rsrcPath + '/close.png'))
        self.closeButton.setIcon(QIcon(rsrcPath + '/close.png'))
        QObject.connect(self.closeButton, SIGNAL("clicked()"), self.close)
        #self.closeButton.clicked.connect(self.close)
        #self.closeButton.connect(self.close())

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        
        self.layout = QGridLayout(self)

        self.layout.addWidget(self.lpix,0,0,1,2, Qt.AlignCenter)       
        self.layout.addWidget(self.label,1,0,1,2)       
        self.layout.addWidget(self.labelblankline,2,0,1,2)       
        #self.layout.addItem(spacerItem,3,0,0,0)

        self.hlayout = QHBoxLayout()
        self.hlayout.addItem(spacerItem)
        self.hlayout.addWidget(self.licenseButton)
        self.hlayout.addWidget(self.closeButton)
        self.layout.addLayout(self.hlayout,3,0,1,2)       

        self.pixmap = QPixmap('images/aboutlogo.png')
        self.lpix.setPixmap(self.pixmap)

        self.setFixedSize(450,501)
        
        self.licenseButton.clicked.connect(self.license)
 
        #QMetaObject.connectSlotsByName(Dialog) 
    def license(self):
        dialog2 = licenseDialog(self)
        dialog2.setAttribute(Qt.WA_DeleteOnClose)
        dialog2.exec_()


class licenseDialog(QDialog):
    def __init__(self,parent=None):
        super(licenseDialog, self).__init__(parent)
        licensetext = """ 
        <HTML>
        <p align="justify" ><font face="Courier"><b>Copyright &copy; 2015 Weverson R. Gomes.</b></p>
        <p><font face="Courier"><b>All rights reserved.</b></p>
        <BR>
        <BR>
        <p align="justify" ><font face="Courier"><b>This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <a href="http://www.gnu.org/licenses">http://www.gnu.org/licenses</a>.</b></font></p>
                      """
 
        self.resize(400, 300)
        self.verticalLayoutWidget = QWidget(self)
        self.verticalLayoutWidget.setGeometry(QRect(4, 10, 391, 281))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout_3 = QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.textBrowser = QTextBrowser(self.verticalLayoutWidget)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.textBrowser.setOpenExternalLinks(True)
        self.verticalLayout_3.addWidget(self.textBrowser)
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pushButton = QPushButton(self.verticalLayoutWidget)
        self.pushButton.setText("&Close")
        self.horizontalLayout.addWidget(self.pushButton)
        self.pushButton.setIcon(QIcon(rsrcPath + '/close.png'))
        QObject.connect(self.pushButton, SIGNAL("clicked()"), self.close)
        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.textBrowser.setText(licensetext)

        #self.retranslateUi(self)
        #QMetaObject.connectSlotsByName(self)

        self.setWindowTitle("License")
