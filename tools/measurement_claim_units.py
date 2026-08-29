"""Frame- and unit-bearing scalar types shared by claim producers."""

from typing import NewType

Degrees = NewType("Degrees", float)
Millimetres = NewType("Millimetres", float)
WorldMillimetres = NewType("WorldMillimetres", float)
