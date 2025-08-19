#!/bin/bash
# Update website, convert Jupyter notebooks to HTML + Markdown, generate index

set -e  # stop on errors

# Directories
NOTEBOOK_DIR="notebooks"
HTML_DIR="notebooks_html"
MD_DIR="_notebooks"

echo "🧹 Cleaning old generated files..."
rm -rf "$HTML_DIR" "$MD_DIR"
mkdir -p "$HTML_DIR" "$MD_DIR"

echo "🔄 Converting notebooks..."

# Convert notebooks
for nb in $NOTEBOOK_DIR/*.ipynb; do
    if [ -f "$nb" ]; then
        base=$(basename "$nb" .ipynb)
        echo "Processing $nb ..."

        # Convert to HTML
        jupyter nbconvert --to html "$nb" --output-dir "$HTML_DIR"

        # Convert to Markdown
        jupyter nbconvert --to markdown "$nb" --output-dir "$MD_DIR"

        # Add YAML front matter if not present
        md_file="$MD_DIR/$base.md"
        if ! grep -q "^---" "$md_file"; then
            sed -i "1i ---\nlayout: page\ntitle: $base\nbasename: $base\n---\n" "$md_file"
        fi
    fi
done

echo "✅ Notebook conversion done."

# Generate notebooks index page
INDEX_FILE="_notebooks/notebooks.md"
echo "🔄 Generating notebooks index page..."
cat > "$INDEX_FILE" <<EOL
---
layout: page
title: Notebooks
permalink: /notebooks/
---

# Jupyter Notebooks

Below is a list of notebooks available in this project.  
Each notebook has three versions:
- **Jekyll page** (rendered nicely in the site)
- **Standalone HTML** (original look, static export)
- **Raw ipynb** (downloadable file)

---

{% for nb in site.notebooks %}
- **{{ nb.title | default: nb.name }}**  
  - [View as page]({{ nb.url }})  
  - [View HTML](../notebooks_html/{{ nb.basename }}.html)  
  - [Download ipynb](../notebooks/{{ nb.basename }}.ipynb)
{% endfor %}
EOL

echo "✅ Notebooks index generated."

# Commit and push changes
git add -A
git commit -m "Update site on $(date)"
git push

echo "🌐 Site updated: https://prasanta13.github.io/compchem_for_newbies/"

