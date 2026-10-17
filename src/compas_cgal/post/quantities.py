"""Units of the emitted program, carried in the type.

Program lengths are exact rationals in the G54 work frame, millimetres (G21).
Feeds and speeds are whole numbers because the writer emits them as integers.
"""

from fractions import Fraction
from typing import NewType

# A length or coordinate in the program's work frame, mm, exact.
ProgramMillimetre = NewType("ProgramMillimetre", Fraction)

# Uniform source-to-program scale, mm per source unit (the gate is unitless), exact.
ScaleFactor = NewType("ScaleFactor", Fraction)

# Feed rate, mm/min. The header sets G94 (feed per minute) explicitly; nothing relies on a controller default.
MillimetrePerMinute = NewType("MillimetrePerMinute", int)

# Spindle speed, rev/min.
RevolutionsPerMinute = NewType("RevolutionsPerMinute", int)

# Digits after the decimal point of every emitted coordinate: the output quantum is 10**-decimals mm.
DecimalPlaces = NewType("DecimalPlaces", int)

# A dwell time, s (G4 P in LinuxCNC), exact.
Seconds = NewType("Seconds", Fraction)
