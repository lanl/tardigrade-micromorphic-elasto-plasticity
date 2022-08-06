# Make bash script more like high-level languages.
set -Eeuxo pipefail
workdir=${PWD}
#============================================================= RUN CPP TESTS ===
# Perform repo tests
cd "build"
ctest --verbose --output-log results.tex
cd ${workdir}
