#!/usr/bin/env bash
# Setup utility to provide the necessary folder structure for data

# find script dir
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# change to script dir, exit if cd fails
cd "$DIR" || exit 1

echo "Starting setup..."

mkdir -p data/images
mkdir -p data/inputs
mkdir -p data/logs
mkdir -p data/outputs/sdf
mkdir -p data/outputs/target_plates

echo "Finished setting up directory structure."

echo "To finish setup:"
echo "- install requirements as specified in requirements.txt"
echo "- place ChemInventory export in data/inputs."
