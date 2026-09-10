# eikonalECGPersonalizedCellZoneRegions

`cellZoneRegions` counterpart of `../eikonalECGPersonalizedNamedRegions`:
identical mesh, stimulus, and `personalizedTemplates` config, but
`ionicHeterogeneity` uses `mode cellZoneRegions` with 3 mesh cell zones
(`system/topoSetDict`, built via `topoSet` in `Allrun`) whose cell-centred
x-ranges reproduce `../eikonalECGPersonalized`'s `transmuralBands`
boundaries (`t=0.3 -> x=0.006`, `t=0.7 -> x=0.014`) as closely as a crisp,
unblended whole-cell partition can.

This is **not** expected to reproduce that tutorial's ECG output
bit-for-bit — `cellZoneRegions` has no blending at all, so cells inside
the `transmuralBands` transition zone get a hard 0/1 assignment here
instead of a smooth weight. What **is** expected to match exactly: each
generated per-region template's peak `Vm`, since
`endocardialCells`/`mCells`/`epicardialCells` baselines resolve to
identical ionic constants regardless of which mode reaches them.

## Stack

Same as `../eikonalECGPersonalizedNamedRegions`, except:

```
ionicHeterogeneity
{
    mode  cellZoneRegions;
    regions
    {
        endocardialCells { cellZone endoZone; }
        mCells           { cellZone midZone;  }
        epicardialCells  { cellZone epiZone;  }
    }
}
```

## Execution

```bash
./Allrun
```

Runs `blockMesh`, then `topoSet` (builds the three cell zones), then
`cardiacFoam`.
