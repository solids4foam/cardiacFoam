"""Standalone regression-equivalence harness for the driverFOAM agent.

Proves the agent reproduces every canonical regression in
``tutorials/Alltest-regression`` by running both the hand-authored path and the
agent-driven path and confirming their numerical outputs agree within each
case's own tolerance. See
``docs/superpowers/plans/2026-07-04-driverfoam-regression-equivalence-harness.md``.
"""
