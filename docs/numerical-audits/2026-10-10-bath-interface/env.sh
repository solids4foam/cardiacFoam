source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
export FOAM_USER_LIBBIN=/private/tmp/bath-mms-audit-20261010/lib
export FOAM_USER_APPBIN=/private/tmp/bath-mms-audit-20261010/bin
export SOLIDS4FOAM_INST_DIR=/private/tmp/bath-mms-audit-20261010/modules/physicsModel
export USE_LIGHTWEIGHT_PHYSICSMODEL=1
export DYLD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$FOAM_LIBBIN:$DYLD_LIBRARY_PATH"
export PATH="$FOAM_USER_APPBIN:$PATH"
export FOAM_LD_LIBRARY_PATH="$DYLD_LIBRARY_PATH"
