# -*- mode: python ; coding: utf-8 -*-
"""
This is the PyInstaller spec file for kkcalc2 GUI application.

It specifies how to build the executable, including data files and filtering out unnecessary DLLs
to reduce the final executable size. Ideally this file is generated using `pyi-makespec` command,
but additional dll filtering is added with python code below.

Note, we also use the `--` additional arguments option to specify the build name.
https://pyinstaller.org/en/stable/spec-files.html#adding-parameters-to-spec-files
"""
import argparse, os
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name") # this is for parameters, use action="store_true" for flags
options = parser.parse_args()

a = Analysis(
    [os.path.join("..", "kkcalc2", "__main__.py")],
    pathex=[],
    binaries=[],
    datas=[
      (os.path.join("..", "kkcalc2", "data"), os.path.join("kkcalc2", "data")),
      (os.path.join("..", "kkcalc2", "asf_database", "ASF.json.gz"), os.path.join("kkcalc2", "asf_database"))],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=2,
)

# -------------------------
## Filter out unnecessary DLLs: # https://github.com/pyinstaller/pyinstaller/issues/4433#issuecomment-671391676
to_keep = []
to_exclude = {'Qt5dbus.dll', 'Qt5Network.dll', 'Qt5Qml.dll', 'Qt5Quick.dll', 'Qt5Svg.dll', 'Qt5WebSockets.dll',
              'libcrypto-3.dll', 'libssl-3.dll', 'sqlite3.dll',
            #   '_sparsetools.cp313-win_amd64.pyd', # Scipy, required
              'cython_special.cp313-win_amd64.pyd', # Scipy
              '_avif.cp313-win_amd64.pyd', # PIL
              'opengl32sw.dll', # PyQt OpenGL Renderer
            #   'libscipy_openblas-48c358d105077551cc9cc3ba79387ed5.dll', # Scipy 32 Bit specific, required
              }

# Iterate through the list of included binaries.
for (dest, source, kind) in a.binaries:
    # Skip anything we don't need.
    if os.path.split(dest)[1] in to_exclude:
        continue
    to_keep.append((dest, source, kind))

# Replace list of data files with filtered one.
a.binaries = to_keep

## Continue spec file as normal.
# -------------------------

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [('O', None, 'OPTION'), ('O', None, 'OPTION')],
    name=options.name if options.name else 'kkcalc2',
    debug=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=[os.path.join("..", "docs", "source", "_static", "logo2.png")]
)
