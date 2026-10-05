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

import json
import tempfile
import webbrowser
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

if TYPE_CHECKING:
    from ...atomic_system import AtomicSystem


def pubchem_search_by_name(name: str) -> list[dict[str, Any]]:
    """Search for CIDs by common name or IUPAC."""
    base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/description/json"
    try:
        response = requests.get(base_url, timeout=10)
        if response.status_code != 200:
            return []

        data = response.json()
        results = []
        for info in data.get("InformationList", {}).get("Information", []):
            cid = info.get("CID")
            if cid:
                # We retrieve the details (Formula, Weight) for display
                details = get_pubchem_details(cid)
                results.append(
                    {
                        "id": str(cid),
                        "common_name": info.get("Title", name),
                        "formula": details.get("formula", "N/A"),
                        "weight": details.get("weight", "N/A"),
                        "smiles": details.get("smiles", "N/A"),
                    }
                )
        return results
    except (requests.RequestException, ValueError, IndexError) as e:
        print(f"PubChem error (Name Search): {e}")
        return []


def pubchem_search_by_formula(formula: str) -> list[dict[str, Any]]:
    """Search for CIDs by chemical formula."""
    base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula}/listkey/json"
    try:
        # Formula search is asynchronous on PubChem (ListKey)
        response = requests.get(base_url, timeout=10)
        if response.status_code != 200:
            return []

        listkey = response.json().get("IdentifierList", {}).get("ListKey")
        # We retrieve the first 20 results for the responsiveness of the UI
        fetch_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/cids/json?maxrecords=20"
        fetch_res = requests.get(fetch_url, timeout=10)
        cids = fetch_res.json().get("IdentifierList", {}).get("CID", [])

        results = []
        for cid in cids:
            details = get_pubchem_details(cid)
            results.append(
                {
                    "id": str(cid),
                    "common_name": details.get("title", "N/A"),
                    "formula": formula,
                    "weight": details.get("weight", "N/A"),
                    "smiles": details.get("smiles", "N/A"),
                }
            )
        return results
    except (requests.RequestException, ValueError, IndexError) as e:
        print(f"PubChem error (Formula Search): {e}")
        return []


def get_pubchem_details(cid: int) -> dict:
    """Retrieves properties. If SMILES is missing, performs a dedicated fallback request."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,Title,IUPACName,CanonicalSMILES/json"

    try:
        r = requests.get(url, timeout=5)
        if r.status_code != 200:
            return {}

        props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        smiles = props.get("ConnectivitySMILES")

        if not smiles:
            fallback_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/json"
            r_s = requests.get(fallback_url, timeout=3)
            if r_s.status_code == 200:
                s_props = r_s.json().get("PropertyTable", {}).get("Properties", [{}])[0]
                smiles = s_props.get("ConnectivitySMILES")

        return {
            "formula": props.get("MolecularFormula"),
            "weight": props.get("MolecularWeight"),
            "smiles": smiles,
            "title": props.get("Title"),
            "iupac": props.get("IUPACName"),
        }
    except (requests.RequestException, ValueError, IndexError) as e:
        print(f"PubChem error (Details): {e}")
        return {}


def get_structure(cid: int, smiles: str = None) -> AtomicSystem | None:
    """Fetch 3D structure from PubChem, fallback to SMILES generation."""
    from cemd.core.atomic_system import AtomicSystem

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    response = requests.get(url, timeout=10)

    if response.status_code == 200:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sdf", delete=False) as f:
            f.write(response.text)
            temp_path = Path(f.name)

        try:
            system = AtomicSystem.from_file(temp_path)
            temp_path.unlink()
            return system
        except Exception as e:
            temp_path.unlink()
            print(f"Error parsing SDF: {e}")

    elif smiles:
        try:
            return AtomicSystem.from_smiles(smiles)
        except Exception as e:
            print(f"SMILES conversion error: {e}")

    return None


def pubchem_sdq_search(query: str, limit: int = 50) -> list[dict]:
    """
    Retrieves exact results from PubChem site via SDQ API.
    """
    # Query cleanup and word separation for "AND"
    words = [w.strip() for w in query.split() if w.strip()]
    if not words:
        return []

    # Construction of the 'where' clause (identical to your URL)
    ands = [{"*": word} for word in words]

    query_params = {
        "download": "*",
        "collection": "compound",
        "order": ["relevancescore,desc"],
        "start": 1,
        "limit": limit,
        "where": {"ands": ands},
    }

    params = {"infmt": "json", "outfmt": "json", "query": json.dumps(query_params)}

    url = "https://pubchem.ncbi.nlm.nih.gov/sdq/sphinxql.cgi"

    try:
        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()  # Throw an error if the query fails

        data = response.json()

        results = []
        for item in data:
            # We map the keys from JSON SDQ to your table format
            results.append(
                {
                    "id": str(item.get("cid", "")),
                    "common_name": item.get("cmpdname", "N/A"),
                    "formula": item.get("mf", "N/A"),
                    "weight": str(item.get("mw", "N/A")),
                    "smiles": item.get("smiles", "N/A"),
                    "iupac": item.get("iupacname", "N/A"),
                }
            )
        return results

    except Exception as e:
        print(f"Error querying SDQ: {e}")
        return []


def explore_pubchem(visible_rows: int = 20) -> Any:
    """Interactive PubChem Explorer aligned with the COD interface."""
    import questionary
    from prompt_toolkit.application import Application
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, ScrollOffsets
    from prompt_toolkit.layout.containers import HSplit, Window
    from prompt_toolkit.layout.controls import FormattedTextControl

    query = questionary.text("Enter name or formula:").ask()
    if not query:
        return None

    print(f"Searching for '{query}' on PubChem...")
    results = pubchem_search_by_name(query)

    if not results:
        print("No results found.")
        return None

    # --------------- Viewer -----------------

    index = 0
    selected = None
    scroll_top = 0

    kb = KeyBindings()

    # Formatting suitable for PubChem columns
    row_format = "{cursor} {id:>10} | {name:<30.30} | {formula:<15.15} | {weight:>10}"

    def render():
        header = row_format.format(
            cursor=" ", id="ID", name="NAME", formula="FORMULA", weight="WEIGHT"
        )
        separator = "-" * len(header)

        lines = [
            f"Found {len(results)} structures  [{index + 1}/{len(results)}]",
            "↑↓ Move   [Enter] Load   [p] PubChem page   [q] Quit",
            "",
            header,
            separator,
        ]

        # Only show the visible window
        visible = results[scroll_top : scroll_top + visible_rows]

        for i, r in enumerate(visible):
            real_index = scroll_top + i
            cursor = "➜" if real_index == index else " "
            lines.append(
                row_format.format(
                    cursor=cursor,
                    id=r["id"],
                    name=r["common_name"],
                    formula=r["formula"],
                    weight=r["weight"],
                )
            )

        # Indicateur de scroll
        if len(results) > visible_rows:
            lines.append(
                f"\n  ↕ {scroll_top + 1}-{min(scroll_top + visible_rows, len(results))} of {len(results)}"
            )

        return "\n".join(lines)

    # Window configuration for scrolling
    text = FormattedTextControl(render)
    window = Window(
        content=text,
        always_hide_cursor=True,
        scroll_offsets=ScrollOffsets(top=1, bottom=1),
        allow_scroll_beyond_bottom=False,
        dont_extend_height=False,
    )

    @kb.add("up")
    def _(e):
        nonlocal index, scroll_top
        if index > 0:
            index -= 1
            if index < scroll_top:
                scroll_top = index
            e.app.invalidate()

    @kb.add("down")
    def _(e):
        nonlocal index, scroll_top
        if index < len(results) - 1:
            index += 1
            if index >= scroll_top + visible_rows:
                scroll_top = index - visible_rows + 1
            e.app.invalidate()

    @kb.add("enter")
    def _(e):
        nonlocal selected
        selected = results[index]
        e.app.exit()

    # Hotkey 'p' to open PubChem page
    @kb.add("p")
    def _(e):
        cid = results[index]["id"]
        webbrowser.open(f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")

    @kb.add("q")
    @kb.add("escape")
    def _(e):
        e.app.exit()

    app = Application(
        layout=Layout(HSplit([window])),
        key_bindings=kb,
        full_screen=True,  # Full screen for more space
    )

    app.run()

    if selected:
        print(f"Loading CID {selected['id']}...")
        return get_structure(int(selected["id"]), selected["smiles"])

    return None
