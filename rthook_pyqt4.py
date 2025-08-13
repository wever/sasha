import sip
import os
import sys

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), "lib"))

sip.setapi(u'QDate', 2)
sip.setapi(u'QDateTime', 2)
sip.setapi(u'QString', 2)
sip.setapi(u'QTextStream', 2)
sip.setapi(u'QTime', 2)
sip.setapi(u'QUrl', 2)
sip.setapi(u'QVariant', 2)
