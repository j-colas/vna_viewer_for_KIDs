# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['/home/jules/Documents/GitRep/vna_viewer_for_KIDs/src/main/python/main.py'],
             pathex=['/home/jules/Documents/GitRep/vna_viewer_for_KIDs/target/PyInstaller'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=['/home/jules/Documents/GitRep/VNA_viewer_for_KIDs/fbsenv/lib/python3.6/site-packages/fbs/freeze/hooks'],
             runtime_hooks=['/tmp/tmpxtwo0ub5/fbs_pyinstaller_hook.py'],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='vna_viewer_4KIDs',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=False,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=False,
               upx_exclude=[],
               name='vna_viewer_4KIDs')
