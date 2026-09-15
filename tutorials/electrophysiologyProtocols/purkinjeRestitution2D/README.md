# purkinjeRestitution2D

A four-branch Purkinje tree coupled to a 50 mm × 50 mm monodomain slab
(Bueno-Orovio epicardial cells, 0.33 mm cells), run over several beats. The
network root is at the left edge; the three terminals are the junctions.

```text
             2  (35, 40) mm
            /
  0 ---- 1 ---- 4  (45, 25) mm
            \
             3  (35, 10) mm
```

The graph in `constant/purkinjeGraph` keeps the branch points as nodes 0-4 and
subdivides each branch into 0.25 mm segments, so the same file serves the
eikonal and the cable network solvers.

## Variants

| Variant | Network solver | Coupler | Protocol |
| --- | --- | --- | --- |
| `antegrade` (default) | `restitutionEikonalSolver1D` | `eikonalMonodomainPvjCoupler`, bidirectional | no stimulus until the network's escape beat at 1.20247 s; root beats at 1.7, 2.1 and 2.4 s |
| `retrograde` | `restitutionEikonalSolver1D` | `eikonalMonodomainPvjCoupler`, bidirectional | tissue stimulus on junction 4 at 0.05 and 0.55 s; no root stimulus |
| `monodomain` | `monodomain1DSolver` | `reactionDiffusionPvjCoupler`, implicit, linear kernel | root current at 0.01 and 0.3 s |

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

## Regression

`regression/regressionTest.sh` runs every variant in parallel and compares
`purkinjeNetwork.dat` against `regression/<variant>.reference`.
For `monodomain` it also runs the graph-only `runPurkinjeGraph` utility.
