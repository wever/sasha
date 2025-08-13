# -*- mode: python -*-

block_cipher = None


a = Analysis(['sasha.py'],
             pathex=['/home/azazel/Downloads/linux_pyqt4/src'],
             binaries=None,
             datas=None,
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=['rthook_pyqt4.py'],
             excludes=None,
             win_no_prefer_redirects=None,
             win_private_assemblies=None,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='sasha',
          debug=False,
          strip=None,
          upx=True,
          console=False , icon='iconsasha.ico')
