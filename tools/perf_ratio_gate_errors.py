"""The ratio gate's refusal family: `GateRefusal`, the base of every named refusal.

A refusal means the gate declined to decide: the registry does not load, the machine is
uncalibrated, the session is not one the registry can judge, a promotion does not qualify.
The CLI (`tools.perf_ratio_gate.main`) maps the whole family to one exit code by type, so a
new refusal joins the class by subclassing this one, never by being listed.

This module sits below every other gate module and imports nothing, so each of them,
`tools.perf_ratio_gate_merge` included, can derive its refusals from it.
"""


class GateRefusal(Exception):
    """The gate declined to decide: the base class of every named ratio gate refusal."""
