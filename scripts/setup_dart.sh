#!/bin/bash
set -e

# Run from the DART root directory
cd $(dirname "$0")/..

echo "Setting up DART build template..."
cp build_templates/mkmf.template.intel.linux build_templates/mkmf.template

echo "Building a low-order model as a compile test..."
cd models/lorenz_63/work

./quickbuild.sh

echo
echo "DART compile test successful."
echo
