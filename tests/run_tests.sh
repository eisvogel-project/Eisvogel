#!/usr/bin/env bash
set -e

BUILD_DIR="$1"
export PHYSICS_TEST_DATA_DIR=${BUILD_DIR}/../tests/physics/data/

for cur in `find ${BUILD_DIR}/tests/ -type f -executable`
do
    echo "---------------------------------------------------------------"
    echo " Running test: '${cur}'"
    echo "---------------------------------------------------------------"
    ${cur}
    echo "---------------------------------------------------------------"
done
