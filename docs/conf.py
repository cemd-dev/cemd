#
# Sphinx configuration for the cemd documentation.
# Full reference: https://www.sphinx-doc.org/en/master/usage/configuration.html

# docs/conf.py

import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(".."))

# ... le reste de votre conf.py ...

# ========================================================================
# Générer les données FF au build
# ========================================================================

# docs/conf.py

import os
import sys

sys.path.insert(0, os.path.abspath(".."))


def generate_ff_data(app):
    """Run scripts to generate force field data from TOML."""
    scripts_dir = Path(__file__).parent / "_scripts"

    # Scripts à exécuter dans l'ordre
    scripts = [
        "db_to_json.py",  # Génère ff_data.json
        "generate_ff_embedded.py",  # Génère ff_viewer.html avec données intégrées
    ]

    for script_name in scripts:
        script_path = scripts_dir / script_name

        if not script_path.exists():
            print(f"⚠️ Script not found: {script_path}")
            continue

        try:
            print(f"🔄 Running {script_name}...")
            result = subprocess.run(
                [sys.executable, str(script_path)],
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent,
                timeout=30,
            )
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(f"⚠️ {result.stderr}")
        except subprocess.TimeoutExpired:
            print(f"⚠️ {script_name} timed out")
        except Exception as e:
            print(f"⚠️ Error running {script_name}: {e}")

    # Copier le JSON dans _build/html/_static/
    src = Path(__file__).parent / "_static" / "ff_data.json"
    dst = Path(__file__).parent / "_build" / "html" / "_static" / "ff_data.json"

    if src.exists():
        try:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            print(f"Copied ff_data.json to {dst}")
        except Exception as e:
            print(f"⚠️ Could not copy ff_data.json: {e}")
    else:
        print(f"⚠️ ff_data.json not found at {src}")

    # Copier le HTML embedded dans _build
    src_html = Path(__file__).parent / "_static" / "ff_viewer.html"
    dst_html = Path(__file__).parent / "_build" / "html" / "_static" / "ff_viewer.html"

    if src_html.exists():
        try:
            dst_html.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_html, dst_html)
            print(f"Copied ff_viewer.html to {dst_html}")
        except Exception as e:
            print(f"⚠️ Could not copy ff_viewer.html: {e}")
    else:
        print(f"⚠️ ff_viewer.html not found at {src_html}")


def setup(app):
    """Sphinx setup hook."""
    app.connect("builder-inited", generate_ff_data)


# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------

project = "cemd"
copyright = "2022-2026, Jérôme Claverie"
author = "Jérôme Claverie"
version = "0.1.0"
release = "0.1.0"

# ---------------------------------------------------------------------------
# General configuration
# ---------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    # Writes a .nojekyll file into the build. Without it GitHub Pages runs
    # the site through Jekyll, which ignores every directory starting with
    # an underscore -- _static and _modules would silently disappear,
    # leaving the site unstyled and the source-code pages missing.
    "sphinx.ext.githubpages",
    "sphinx_design",
    "numpydoc",
]

templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "sphinx"
todo_include_todos = True
add_module_names = False

# ---------------------------------------------------------------------------
# Autodoc / Autosummary
# ---------------------------------------------------------------------------

autodoc_typehints = "none"
autodoc_member_order = "bysource"

# Heavy optional dependencies — mock them so Sphinx does not need to import them
autodoc_mock_imports = [
    "pymatgen",
    "rdkit",
]

autodoc_default_options = {
    "show-inheritance": True,
    "undoc-members": False,
    # Exclude instance attributes already documented in the class docstring
    # to avoid them appearing twice in the rendered page.
    "exclude-members": (
        "atoms, bonds, angles, dihedrals, impropers, velocities"
    ),
}

autosummary_generate = True
autosummary_imported_members = True

# ---------------------------------------------------------------------------
# NumpyDoc
# ---------------------------------------------------------------------------

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_attributes_as_param_list = False

# ---------------------------------------------------------------------------
# Intersphinx — link to external docs
# ---------------------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "MDAnalysis": ("https://docs.mdanalysis.org/stable", None),
}


# ---------------------------------------------------------------------------
# HTML output — PyData Sphinx Theme
# ---------------------------------------------------------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "logo": {
        "image_light": "_static/images/logo/cemd_logo.svg",
        "image_dark": "_static/images/logo/cemd_logo_dark.svg",
    },
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/cemd-dev/cemd",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}

html_show_sourcelink = False

html_static_path = ["_static"]

html_css_files = [
    "https://fonts.googleapis.com/icon?family=Material+Icons",
    "custom.css",
]

html_sidebars = {
    "api/index": [],
    "api/atomic_system": [],
    "user_guide/index": [],
    "user_guide/build_guide": [],
    "user_guide/analysis_guide": [],
    "user_guide/ff_database": [],
    "getting_started": [],
    "installation": [],
}

# ---------------------------------------------------------------------------
# LaTeX output
# ---------------------------------------------------------------------------

# latex_documents = [
#     (master_doc, 'cemd.tex', 'cemd Documentation', author, 'manual'),
# ]


# todo_include_todos = True

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "cemdsdoc"

# ---------------------------------------------------------------------------
# Other build (HTMLHelp, man, Texinfo, epub)
# ---------------------------------------------------------------------------

htmlhelp_basename = "cemddoc"

man_pages = [
    (master_doc, "cemd", "cemd Documentation", [author], 1),
]

texinfo_documents = [
    (
        master_doc,
        "cemd",
        "cemd Documentation",
        author,
        "cemd",
        "Atomistic model builder for LAMMPS.",
        "Miscellaneous",
    ),
]

epub_title = project
epub_exclude_files = ["search.html"]
