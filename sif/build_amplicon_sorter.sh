#!/bin/bash

# Build script for amplicon_sorter Singularity container
# Usage: ./sif/build_amplicon_sorter.sh

set -e

echo "Building Singularity container for amplicon_sorter..."
singularity build --fakeroot sif/amplicon_sorter_v2.sif sif/amplicon_sorter.def

echo "Build completed successfully. Container: sif/amplicon_sorter_v2.sif"