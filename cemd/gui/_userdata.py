"""Where the interface keeps the files it writes.

The package directory is not one of them. It works in a source tree and
fails everywhere else: an installed package lives under ``site-packages``,
which is shared between the users of a machine, replaced wholesale on
every upgrade, and read-only as often as not. Settings written there are
silently lost at best and refused at worst.

Qt already knows the per-user location each platform expects, so ask it
rather than guessing. The generic locations are used and the application
directory appended here, because the application-specific ones depend on
``QCoreApplication.applicationName()`` being set -- which it is not when a
dialog is built outside :func:`cemd.gui.main_window.main`, as the tests do.
"""

from __future__ import annotations

from pathlib import Path

from PySide6 import QtCore

APP_DIR = "cemd"


def _resolve(location: QtCore.QStandardPaths.StandardLocation, fallback: str) -> Path:
    """Return the per-user *location*, falling back if Qt reports none."""
    root = QtCore.QStandardPaths.writableLocation(location)
    return Path(root or Path.home() / fallback) / APP_DIR


def config_dir() -> Path:
    """Directory holding the user's saved preferences."""
    return _resolve(QtCore.QStandardPaths.GenericConfigLocation, ".config")


def cache_dir() -> Path:
    """Directory holding what can be thrown away and fetched again."""
    return _resolve(QtCore.QStandardPaths.GenericCacheLocation, ".cache")


def config_file(name: str = "config.json") -> Path:
    """Path of a preferences file, its directory created if needed."""
    directory = config_dir()
    directory.mkdir(parents=True, exist_ok=True)
    return directory / name


def cache_file(name: str) -> Path:
    """Path of a cache file, its directory created if needed."""
    directory = cache_dir()
    directory.mkdir(parents=True, exist_ok=True)
    return directory / name


def legacy_config_file() -> Path:
    """The place earlier versions wrote to, inside the package itself.

    Read once, so that upgrading does not discard preferences that are
    already there.
    """
    return Path(__file__).resolve().parent / "config.json"
