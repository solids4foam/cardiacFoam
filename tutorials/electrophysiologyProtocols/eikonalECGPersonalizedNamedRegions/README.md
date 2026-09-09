# eikonalECGPersonalizedNamedRegions

`namedRegions` counterpart of `../eikonalECGPersonalized`: identical mesh,
stimulus, and `personalizedTemplates` config, but `ionicHeterogeneity` uses
`mode namedRegions` with 3 regions whose ranges reproduce
`eikonalECGPersonalized`'s `transmuralBands` boundaries (`0/0.3/0.7/1`)
exactly. Its `postProcessing/eikonalECG.dat` is expected to match that
tutorial's output to numerical precision — an end-to-end version of the
N=3 `transmuralBands`/`namedRegions` equivalence proved at the unit level
by `src/ionicModels/tests/test_ionic_heterogeneity_synthesize_transmural_regions.py`.

## Stack

Same as `../eikonalECGPersonalized`, except:

```
ionicHeterogeneity
{
    mode  namedRegions;
    regions
    {
        endocardialCells { range (0   0.3); }
        mCells           { range (0.3 0.7); }
        epicardialCells  { range (0.7 1  ); }
    }
}
```

## Execution

```bash
./Allrun
```
