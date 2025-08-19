#!/bin/bash

#========================================
# update_site.sh
# Converts notebooks to HTML/Markdown
# and pushes to GitHub Pages
#========================================

set -e  # Exit immediately if a command fails

# Directories
NOTEBOOKS_DIR="notebooks"
HTML_DIR="${NOTEBOOKS_DIR}_html"
MD_DIR="_notebooks"

# Create output directories if they don't exist
mkdir -p "$HTML_DIR"
mkdir -p "$MD_DIR"

# Clean old generated files
echo "🧹 Cleaning old generated files..."
rm -f "$HTML_DIR"/*.html
rm -f "$MD_DIR"/*.md

# Convert notebooks
echo "🔄 Converting notebooks..."
for nb in "$NOTEBOOKS_DIR"/*.ipynb; do
    echo "Processing $nb ..."

    base=$(basename "$nb" .ipynb)

    # Convert to HTML
    jupyter nbconvert \
        --to html \
        --output "$base" \
        --output-dir "$HTML_DIR" \
        "$nb"

    # Convert to Markdown
    jupyter nbconvert \
        --to markdown \
        --output "$base" \
        --output-dir "$MD_DIR" \
        "$nb"
done

echo "✅ Conversion completed."

# Git auto-commit and push
echo "📦 Adding changes to Git..."
git add -A

commit_msg="Update site: $(date +'%Y-%m-%d %H:%M:%S')"
echo "📝 Commit message: $commit_msg"
git commit -m "$commit_msg" || echo "⚠️ Nothing to commit."

echo "🚀 Pushing to GitHub..."
git push origin main

echo "🎉 Site updated successfully!"

