#!/bin/bash
set -e

# Build the Jupyter Book
cd mybook
echo "📚 Building Jupyter Book..."
jupyter-book build .

# Ensure GitHub Pages won't run Jekyll on the output
touch _build/html/.nojekyll
cd ..

# Prepare gh-pages branch
echo "🔀 Switching to gh-pages..."
if git show-ref --verify --quiet refs/heads/gh-pages; then
  git checkout gh-pages
else
  git checkout --orphan gh-pages
  # Remove all files from index (if any)
  git rm -rf . >/dev/null 2>&1 || true
fi

# Wipe old site contents
find . -maxdepth 1 ! -name '.' ! -name '..' ! -name '.git' -exec rm -rf {} +

# Copy fresh site
cp -r mybook/_build/html/* .

# Ask for commit message
read -rp "Enter commit message for site update: " MSG
git add .
git commit -m "${MSG:-Update site}"
git push -u origin gh-pages

# Back to main branch
git checkout main
echo "✅ Deployed to gh-pages. Visit: https://prasanta13.github.io/compchem_for_newbies/"
