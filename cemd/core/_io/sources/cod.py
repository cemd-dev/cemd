#
# This file is part of the CEMD distribution
# Copyright (c) 2022-2026 Jérôme Claverie.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import annotations

import re
import warnings
from io import StringIO
from typing import TYPE_CHECKING, Any

import requests
from pymatgen.io.cif import CifParser

if TYPE_CHECKING:
    from ...atomic_system import AtomicSystem


def formula_to_elements(formula: str) -> list[dict[str, Any]]:
    """
    Transforms a chemical formula into a list of unique elements.
    Example: "Ca(OH)2" -> ['Ca', 'O', 'H']
              "Al2O3" -> ['Al', 'O']
    """

    return list(dict.fromkeys(re.findall(r"[A-Z][a-z]?", formula)))


def cod_search_by_elements(elements_list: list) -> list[dict[str, Any]]:
    """
    EXCLUSIVE search: only finds compound structures
    STRICTLY items provided.
    """
    # Cleaning and sorting
    elements = sorted([el.strip() for el in elements_list if el.strip()])
    num_el = len(elements)

    base_url = "https://www.crystallography.net/cod/result.php"
    params = {
        "el[]": elements,
        "strictmin": num_el,  # Minimum n elements
        "strictmax": num_el,  # Maximum n elements
        "format": "json",
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code != 200:
            return []

        results = response.json()
        if not isinstance(results, list):
            return []

        cleaned_results = []
        # Define a set of our elements in capital letters for comparison
        target_set = set(el.upper() for el in elements)

        for r in results:
            formula = r.get("formula", "")
            if not formula:
                continue

            found_elements = set(formula_to_elements(formula))
            found_set_upper = set(el.upper() for el in found_elements)

            if found_set_upper == target_set:
                cleaned_results.append(
                    {
                        "id": r.get("id") or r.get("file"),
                        "formula": formula,
                        "spacegroup": r.get("sg", "N/A"),
                        "a": float(r.get("a", 0)),
                        "b": float(r.get("b", 0)),
                        "c": float(r.get("c", 0)),
                        "alpha": float(r.get("alpha", 0)),
                        "beta": float(r.get("beta", 0)),
                        "gamma": float(r.get("gamma", 0)),
                        "common_name": r.get("mineral") or r.get("compound") or "N/A",
                        "doi": r.get("doi", "N/A"),
                    }
                )
        return cleaned_results

    except Exception as e:
        print(f"Search error: {e}")
        return []


def cod_search_by_name(mineral_name: str) -> list[dict[str, Any]]:
    """
    Search the COD by mineral name or free text.
    Returns a formatted list with individual lattice parameters.
    """
    base_url = "https://www.crystallography.net/cod/result.php"
    params = {"text": mineral_name, "format": "json"}

    try:
        response = requests.get(base_url, params=params, timeout=30)

        if response.status_code == 200:
            results = response.json()
            if not isinstance(results, list):
                return []

            cleaned_results = []
            for r in results:
                # extract the values with a secure conversion to float
                # If the value does not exist or is not a number, we set 0.0
                def to_float(val):
                    try:
                        return float(val)
                    except (TypeError, ValueError):
                        return 0.0

                cleaned_results.append(
                    {
                        "id": r.get("id") or r.get("file"),
                        "formula": r.get("formula", "N/A"),
                        "spacegroup": r.get("sg", "N/A"),
                        # Separate storage to allow dynamic formatting
                        "a": to_float(r.get("a")),
                        "b": to_float(r.get("b")),
                        "c": to_float(r.get("c")),
                        "alpha": to_float(r.get("alpha")),
                        "beta": to_float(r.get("beta")),
                        "gamma": to_float(r.get("gamma")),
                        "common_name": r.get("mineral") or r.get("compound") or "N/A",
                        "doi": r.get("doi", "N/A"),
                    }
                )
            return cleaned_results

    except Exception as e:
        print(f"Connection error (Name Search): {e}")

    return []


def cod_search_by_id(cod_id: int) -> list[dict[str, Any]]:
    """Performs a direct search on a specific ID."""
    base_url = "https://www.crystallography.net/cod/result.php"
    params = {"id": cod_id, "format": "json"}
    try:
        r = requests.get(base_url, params=params, timeout=10)
        if r.status_code == 200:
            data = r.json()
            # Format to match 'cleaned_results' format
            return [
                {
                    "id": d.get("id"),
                    "formula": d.get("formula", "N/A"),
                    "spacegroup": d.get("sg", "N/A"),
                    "a": float(d.get("a", 0)),
                    "b": float(d.get("b", 0)),
                    "c": float(d.get("c", 0)),
                    "alpha": float(d.get("alpha", 0)),
                    "beta": float(d.get("beta", 0)),
                    "gamma": float(d.get("gamma", 0)),
                    "common_name": d.get("mineral") or d.get("compound") or "N/A",
                    "doi": d.get("doi", "N/A"),
                }
                for d in data
            ]
    except (requests.RequestException, ValueError) as e:
        print(f"Connection error (ID Search): {e}")
    return []


def get_structure_by_cod_id(cod_id: int) -> AtomicSystem:
    """
    Retrieves a structure without the need for MySQL.
    We download the raw CIF and parse it with Pymatgen.
    """

    from ...atomic_system import AtomicSystem

    # Direct URL of CIF file
    url = f"https://www.crystallography.net/cod/{cod_id}.cif"
    # try:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    parser = CifParser(StringIO(response.text))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        pmg_structure = parser.parse_structures()[0]
    return AtomicSystem.from_pmg(pmg_structure)


def explore_cod(visibles_rows: int = 20) -> AtomicSystem | None:
    """Interactive COD explorer."""

    import webbrowser

    import questionary
    from prompt_toolkit.application import Application
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, ScrollOffsets
    from prompt_toolkit.layout.containers import HSplit, Window
    from prompt_toolkit.layout.controls import FormattedTextControl

    mode = questionary.select(
        "Search by:",
        choices=[
            "Mineral Name",
            "Chemical Formula",
            "COD ID",
        ],
    ).ask()

    if mode is None:
        return None

    query = questionary.text("Enter your query:").ask()

    if not query:
        return None

    # ----------------Search ----------------

    if mode == "Mineral Name":
        results = cod_search_by_name(query)

    elif mode == "Chemical Formula":
        results = cod_search_by_elements(formula_to_elements(query))

    else:
        return get_structure_by_cod_id(query)

    if not results:
        print("No structures found.")
        return None

    # ---------------Viewer -----------------

    index = 0
    selected = None
    scroll_top = 0

    kb = KeyBindings()

    row_format = "{cursor} {id:>8} | {name:<20.20} | {formula:<18.18} | {lattice:<37.37} | {doi:^3}"

    def render():
        header = row_format.format(
            cursor=" ",
            id="ID",
            name="NAME",
            formula="FORMULA",
            lattice="LATTICE PARAMETERS",
            doi="DOI",
        )
        separator = "-" * len(header)

        lines = [
            f"Found {len(results)} structures  [{index + 1}/{len(results)}]",
            "↑↓ Move   [Enter] Load   [d] DOI   [c] COD page   [q] Quit",
            "",
            header,
            separator,
        ]

        # Show only the visible window
        visible = results[scroll_top : scroll_top + visibles_rows]

        for i, r in enumerate(visible):
            real_index = scroll_top + i
            cursor = "➜" if real_index == index else " "
            doi_val = "✓" if r.get("doi") and r.get("doi") != "N/A" else "-"
            lat = f"{r.get('a', 0):>5.1f} {r.get('b', 0):>5.1f} {r.get('c', 0):>5.1f}"
            ang = f"{r.get('alpha', 0):>3.0f}° {r.get('beta', 0):>3.0f}° {r.get('gamma', 0):>3.0f}°"

            lines.append(
                row_format.format(
                    cursor=cursor,
                    id=str(r.get("id", "N/A")),
                    name=str(r.get("common_name") or "N/A"),
                    formula=str(r.get("formula", "N/A")),
                    lattice=f"{lat} | {ang}",
                    doi=doi_val,
                )
            )

        # Scroll indicator
        if len(results) > visibles_rows:
            lines.append(
                f"\n  ↕ {scroll_top + 1}-{min(scroll_top + visibles_rows, len(results))} of {len(results)}"
            )

        return "\n".join(lines)

    # Window Definition
    text = FormattedTextControl(render)
    window = Window(
        content=text,
        always_hide_cursor=True,
        scroll_offsets=ScrollOffsets(top=1, bottom=1),
        allow_scroll_beyond_bottom=False,
        dont_extend_height=False,
    )

    @kb.add("up")
    def _(event):
        nonlocal index, scroll_top
        if index > 0:
            index -= 1
            if index < scroll_top:
                scroll_top = index
            event.app.invalidate()

    @kb.add("down")
    def _(event):
        nonlocal index, scroll_top
        if index < len(results) - 1:
            index += 1
            if index >= scroll_top + visibles_rows:
                scroll_top = index - visibles_rows + 1
            event.app.invalidate()

    @kb.add("enter")
    def _(event):
        nonlocal selected
        selected = results[index]
        event.app.exit()

    @kb.add("d")
    def _(event):
        doi = results[index].get("doi")
        if doi and doi != "N/A":
            webbrowser.open(f"https://doi.org/{doi}")

    @kb.add("c")
    def _(event):
        cod = results[index]["id"]
        webbrowser.open(f"https://www.crystallography.net/cod/{cod}.html")

    @kb.add("q")
    @kb.add("escape")
    def _(event):
        event.app.exit()

    app = Application(
        layout=Layout(HSplit([window])),
        key_bindings=kb,
        full_screen=True,
    )

    app.run()

    if selected is None:
        return None

    print(f"Loading structure {selected['id']}...")

    return get_structure_by_cod_id(selected["id"])
