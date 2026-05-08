#!/bin/bash

# Find the correctionlib library path
export CORR_BASE=$(python3 -c "import correctionlib; print(correctionlib.__path__[0])")

# Compile
g++ -shared -fPIC -o libJECUtils.so src/JECUtils.cc \
    $(root-config --cflags --libs) \
    -D_GLIBCXX_USE_CXX11_ABI=0 \
    -I${CORR_BASE}/include -Iinterface \
    -L${CORR_BASE}/lib -lcorrectionlib \
    -Wl,-rpath,${CORR_BASE}/lib

echo "Libraries compiled! You are ready to go."

