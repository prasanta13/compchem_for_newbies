#!/bin/bash
# Update website and convert notebooks to HTML

set -e  # stop on errors

# Directory for HTML notebooks
NOTEBOOK_HTML_DIR="notebooks_html"
mkdir -p $NOTEBOOK_HTML_DIR

# Convert all notebooks in notebooks/ to HTML
for nb in notebooks/*.ipynb; do
    if [ -f "$nb" ]; then
        jupyter nbconvert --to html "$nb" --output-dir "$NOTEBOOK_HTML_DIR"
    fi
done

# Commit changes
git add -A
git commit -m "Update site on $(date)"
##git push
git push -u origin main

echo "✅ Site updated: https://prasanta13.github.io/compchem_for_newbies/"

