#!/bin/bash
set -e

# Build the Jupyter Book
cd mybook || exit
jupyter-book build .

# Go back
cd ..

# Ask user for commit message
read -rp "Enter commit message: " msg

# Stage and commit everything
git add .
git commit -m "$msg"
git push origin main

echo "✅ Book built and pushed. Pages will update at:"
echo "   https://prasanta13.github.io/compchem_for_newbies/"

