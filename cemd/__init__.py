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

from importlib.metadata import PackageNotFoundError, version

from .core.atomic_system import AtomicSystem

try:
    __version__ = version("cemd")
except PackageNotFoundError:  # running from a source tree, never installed
    __version__ = "unknown"

__author__ = "Jérôme Claverie"

__copyright__ = "Copyright (c) 2022-2026 Jérôme Claverie"

__license__ = "GPL-3.0"

__all__ = ["AtomicSystem", "__version__"]
