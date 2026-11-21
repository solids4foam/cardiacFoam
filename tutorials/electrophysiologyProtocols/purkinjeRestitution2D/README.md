# purkinjeRestitution2D

## What it does

A four-branch Purkinje tree coupled to a 50 mm × 50 mm monodomain slab
(Bueno-Orovio epicardial cells, 0.33 mm cells), run over several beats in
three variants: `antegrade` and `retrograde` drive the network with
`restitutionEikonalSolver1D`, `monodomain` with `monodomain1DSolver`. The
network root is at the left edge; the three terminals are the junctions.
Each run writes the network's node values to
`postProcessing/purkinjeNetwork.dat`.

- **Exercises:** the compiled restitution defaults of
  `restitutionEikonalSolver1D` over several beats (escape beat, slower
  conduction at a shorter diastolic interval, block at the root);
  `eikonalMonodomainPvjCoupler` in bidirectional mode, with activation
  entering the network from the root (`antegrade`) and from the tissue
  (`retrograde`); `reactionDiffusionPvjCoupler` with the implicit coupling
  scheme and the linear kernel (`monodomain`); and the graph-only
  `runPurkinjeGraph` utility. Only this case checks these against reference
  values.
- **Use it for:** a starting point for Purkinje-myocardium coupling on a
  slab; checking capture, conduction slowing and block of a restitution
  network at chosen pacing intervals; comparing eikonal and cable networks on
  the same graph.

## Regression

Tier: `full` (marker `regression/full-only`).

`regression/regressionTest.sh` runs `Allrun <variant> parallel` (4
subdomains, `scotch`) for `retrograde`, `monodomain` and `antegrade`, and
compares `postProcessing/purkinjeNetwork.dat` with frozen values from earlier
runs in `regression/<variant>.reference`: node activation times for
`retrograde` (8 rows) and `antegrade` (12 rows), node `Vm` and junction
coupling current for `monodomain` (8 rows). `retrograde` stops at 0.6 s
instead of 0.9 s, once the second beat has crossed the network; `antegrade`
runs on a 75 x 75 slab instead of 150 x 150, since its checks are network
activation times, which the slab resolution does not change. After the
`monodomain` run it also runs `runPurkinjeGraph` on the network alone and
checks 4 rows of `regression/monodomain.graphUtility.reference`: node `Vm`
matches the coupled run and the junction current is zero. Tolerances are
1e-6 for activation times and `Vm`, and 5 A/m³ for the coupled junction
current. `--check-only` is not supported, since each variant needs its own
run. Run it on its own from the case root:

```bash
./regression/regressionTest.sh
```

## Variants

| Variant | Network solver | Coupler | Protocol |
| --- | --- | --- | --- |
| `antegrade` (default) | `restitutionEikonalSolver1D` | `eikonalMonodomainPvjCoupler`, bidirectional | no stimulus until the network's escape beat at 1.20247 s; root beats at 1.7, 2.1 and 2.4 s |
| `retrograde` | `restitutionEikonalSolver1D` | `eikonalMonodomainPvjCoupler`, bidirectional | tissue stimulus on junction 4 at 0.05 and 0.55 s; no root stimulus |
| `monodomain` | `monodomain1DSolver` | `reactionDiffusionPvjCoupler`, unidirectional (network to tissue), implicit, linear kernel | root current at 0.01 and 0.3 s |

`antegrade` exercises the restitution solver's compiled defaults. After the
escape beat, the 1.7 s beat arrives at an interval of 0.4975 s and the 2.1 s
beat at 0.4 s. Both are captured, and the second conducts more slowly because
its diastolic interval is shorter. The 2.4 s beat falls below
`apdNominal + minimumDI90` (0.353 s) and is blocked at the root.

`retrograde` has no network stimulus. The tissue activates junction 4, the
wave climbs the network to the root and returns to the tissue through
junctions 2 and 3. The second tissue beat is captured at a shorter diastolic
interval and travels up the tree more slowly than the first.

`monodomain` integrates the cable on every graph edge and exports node `Vm`
and the junction coupling current.

The restitution constants and how they were measured are described in
`../cableProtocol/monodomain1DCableCV/Purkinje_S1_S2_Calibration.md`.

## Running

```bash
./Allrun                      # antegrade
./Allrun retrograde
./Allrun monodomain parallel  # 4 subdomains, system/decomposeParDict
./Allclean
```

`Allrun` links `constant/electroProperties` and `system/controlDict` to the
chosen variant.

## Outputs

`postProcessing/purkinjeNetwork.dat` holds one row per write time. The
restitution variants write each node's latest activation time. The
`monodomain` variant writes each junction's coupling current, then each
node's `Vm`. `postProcessing/purkinjeNetworkVTK/` holds the network as VTK.

## Notes on dictionary entries

**Mesh.** The 50 x 50 x 0.1 mm slab has 150 x 150 x 1 cells, deltaX = 0.33 mm
(the same mesh as `rotorInstability`).

**Purkinje graph.** `constant/purkinjeGraph` holds four straight branches from
the root, subdivided into 0.25 mm segments so `monodomain1DSolver` has a cable
to integrate; the same file serves the eikonal and the cable network solvers. Nodes 0-4 are the branch points and keep their indices; nodes 5
onwards are interior. The file gives coordinates in metres; the diagram gives
them in mm.

```
             2  (35, 40)
            /
  0 ---- 1 ---- 4  (45, 25)
            \
             3  (35, 10)
```

`conductionEdges` entries are `(nodeA nodeB length conductance)`.

**Variants.**

- `electroProperties.antegrade`: no tissue stimulus; every myocardial
  activation enters through a junction. `apdNominal`, `minimumDI90` and
  `escapeInterval` use the compiled defaults. The escape occurs at 1.20247 s;
  root beats at intervals of 0.4975 s and 0.4 s are captured, and the beat at
  0.3 s falls below `apdNominal + minimumDI90`.
- `electroProperties.retrograde`: two tissue beats 0.5 s apart, each on the
  junction at node 4.
- `electroProperties.monodomain`: current is injected at the root node over
  each stimulus window.
