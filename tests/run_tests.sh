#!/usr/bin/env bash
set -e

BUILD_DIR="$1"

for cur in `find ${BUILD_DIR}/tests/ -type f -executable`
do
    echo "---------------------------------------------------------------"
    echo " Running test: '${cur}'"
    echo "---------------------------------------------------------------"
    ./${cur}
    echo "---------------------------------------------------------------"
done
