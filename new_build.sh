# USAGE:
#
# ./new_build.sh cmake_build_type

# Make bash script more like high-level languages.
set -Eeuxo pipefail

# Get this scripts file name
script=`basename "$0"`

# Parse arguments
if [ "$#" -ne 1 ]; then
    echo "${script} USAGE:"
    echo "./${script} cmake_build_type"
    echo "    cmake_build_type: string for the CMake config -DCMAKE_BUILD_TYPE=<string> option"
    exit 1
fi
cmake_build_type=$1

# Debugging
whoami
ls -l $HOME/include || true
ls -l $HOME/.local/include || true

# Find cmake3 executable
if [ -x "$(command -v cmake3)" ]; then
    cmake_exec=$(command -v cmake3)
elif [ -x "$(command -v cmake)" ]; then
    cmake_exec=$(command -v cmake)
else
    echo "Could not find cmake executable"
    exit 3
fi

# Clean and build repo tests
rm -rf build/
mkdir build
cd build
${cmake_exec} .. -DCMAKE_BUILD_TYPE=${cmake_build_type}
${cmake_exec} --build . --verbose
