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
