# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'license.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
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

class licenseDialog(object):
    def setupUi(self, Dialog2):
        Dialog2.setObjectName(_fromUtf8("Dialog2"))
        Dialog2.resize(400, 300)
        self.verticalLayoutWidget = QtGui.QWidget(Dialog2)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(4, 10, 391, 281))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.textBrowser = QtGui.QTextBrowser(self.verticalLayoutWidget)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.verticalLayout_3.addWidget(self.textBrowser)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pushButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.horizontalLayout.addWidget(self.pushButton)
        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog2)
        QtCore.QMetaObject.connectSlotsByName(Dialog2)

    def retranslateUi(self, Dialog2):
        Dialog2.setWindowTitle(_translate("Dialog2", "Dialog2", None))
        self.pushButton.setText(_translate("Dialog2", "&Close", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

