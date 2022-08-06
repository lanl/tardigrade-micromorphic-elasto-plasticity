# Make bash script more like high-level languages.
set -Eeuxo pipefail
# Find cmake3 executable
if [ -x "$(command -v cmake3)" ]; then
    cmake_exec=$(command -v cmake3)
elif [ -x "$(command -v cmake)" ]; then
    cmake_exec=$(command -v cmake)
else
    echo "Could not find cmake executable"
    exit 3
fi
# Get current conda environment information or exit on error
conda_env_path=$(conda info | grep "active env location" | cut -f 2 -d :)
# Change to build directory and run cmake install
cd "build"
${cmake_exec} --install . --prefix ${conda_env_path}
