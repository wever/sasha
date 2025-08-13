# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aboutDialog2.ui'
#
# Created by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4.QtGui import  QFileDialog,QFileDialog,QMessageBox,QApplication,QMainWindow,QPrinter,QAbstractPrintDialog,QSplashScreen,QMovie,QPixmap,QDialog#,QDesktopWidget,QWidget,QFont,QDialog,QIcon,QVBoxLayout,QHBoxLayout,QTabWidget,QGridLayout,QLabel,QSizePolicy,QSpacerItem,QComboBox,QSpinBox,QCheckBox,QTextEdit,QPushButton,QMenuBar,QMenu,QStatusBar,QAction,QKeySequence,qApp
 
from PyQt4.QtCore import QObject,QThread,SIGNAL,pyqtSignal,QCoreApplication,Qt,QRect,QFileInfo,QFile,QTextStream,QSettings,QMetaObject

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(400, 301)
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(10, 100, 381, 151))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayoutWidget = QtGui.QWidget(Dialog)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(0, 260, 391, 41))
        self.horizontalLayoutWidget.setObjectName(_fromUtf8("horizontalLayoutWidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pushButton = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.horizontalLayout.addWidget(self.pushButton)
        self.pushButton_2 = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.horizontalLayout.addWidget(self.pushButton_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "dfgtedfgsdgsdfgsdfgsd fgdsfgdsf", None))
        self.pushButton.setText(_translate("Dialog", "&License", None))
        self.pushButton_2.setText(_translate("Dialog", "&Close", None))

