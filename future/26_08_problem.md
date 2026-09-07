# driverFOAM dictionary scanner question — `cm`, `chi`, and `stale_paths` (2026-08-26)

## Question

Why do some catalogue paths appear under `stale_paths`, while keys such as
`cm` and `chi` do not, even though the myocardium code reads them through the
same kind of OpenFOAM dictionary object? Does this indicate that strict plan
or sweep validation is unreliable?

## Verified conclusion

The strict planning and sweep workflow is currently valid. The issue is a
limitation in the source-audit scanner's precision, not a demonstrated solver
or catalogue failure.

The scanner does not reconstruct the full dictionary path from C++ source. It
recognises literal method calls such as:

```cpp
dict.lookup("key");
dict.get<scalar>("key");
dict.subDict("block");
```

It then compares bare names globally. It does not track which dictionary a
receiver represents, and it does not currently understand constructor-based
or variable-based reads.

## The `cm`/`chi` example

The myocardium domain reads its physical coefficients through OpenFOAM's
dimensioned-value constructor:

```cpp
chi_("chi", ..., electroProperties_)
Cm_("cm", ..., electroProperties_)
```

See `src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C`.
These reads are not visible to the scanner.

The conduction-system domain separately contains literal reads:

```cpp
chi_(coeffsDict_.get<scalar>("chi"))
Cm_(coeffsDict_.get<scalar>("cm"))
```

See `src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C`.
Those reads make the bare names `chi` and `cm` visible globally, so the
scanner treats both catalogue paths as observed. It has not proved that the
myocardium paths were observed specifically.

This is a source-audit precision limitation. It does not make either runtime
read incorrect: the myocardium constructor receives `electroProperties_`,
while the conduction model receives `purkinjeGraphModelCoeffs`.

## What the `stale_paths` mean

`stale_paths` means that a catalogue leaf was not observed by this approximate
scanner anywhere in the source tree. It does not mean that the runtime key is
unused or invalid.

The current examples fall into three groups:

- `solver`, `absTol`, `relTol`, and `maxSteps` are passed through the ionic
  model to OpenFOAM's `ODESolver`, which reads them upstream.
- `c0` is read with `dimensionedScalar("c0", electroProperties)`, so the key
  is supplied as a constructor argument rather than a literal dictionary
  method call.
- conductivity and stimulus values use shared helpers. Conductivity names
  are held in `conductivityFieldSpec.dictionaryEntry`, while stimulus
  duration and intensity use `dimensionedScalar` constructors.

These are legitimate dictionary access patterns, but the scanner does not
parse them.

## Why strict plan and sweep still pass

The strict planner uses the catalogue's full paths when generating and
checking case dictionaries. The case-key checker is path-aware and matches
wildcard domain names by position. Its unmatched-key diagnostics are
intentionally warnings, because OpenFOAM and external libraries may own keys
outside the plugin catalogue.

The C++/catalogue scanner is a separate aggregate vocabulary check. It is
allowlist-backed and does not claim exact source-to-path attribution.

Verified in the current tree:

- strict dictionary scanner: `status: ok`;
- dictionary scanner, strict-planning, case-dictionary, sweep-plan, and
  runtime-enum tests: `39 passed, 14 subtests passed`;
- catalogue, builder, and validation tests: `167 passed`;
- plugin architecture tests: `3 passed`;
- OpenFOAM v2412 lightweight native library build completed successfully.

## Decision

Use two complementary sources of truth, with different responsibilities:

- The `DictEntry` catalogue is authoritative for the driver-facing dictionary
  contract: paths, exposed types, constraints, applicability, and values the
  builder may emit or override.
- Each solver and shared C++ helper is authoritative for runtime behaviour:
  accepted reads, fallback defaults, units, and the effect of omitting a key.
  A catalogue `typical_value` is a construction default/template and must not
  be assumed to be the same as the solver's C++ fallback.

`source_refs` are provenance links between those two layers. They are useful
for traceability and file-existence checks, but they do not drive dictionary
construction. The C++ scanner is an optional maintenance guard; it should
detect possible drift without becoming a third source of truth. Its allowlist
records reviewed scanner limitations, including constructor reads, helper
indirection, and upstream OpenFOAM pass-throughs.

Do not delete the allowlist or redesign the scanner as part of this issue. The
only allowlist cleanup needed here was removal of the unused
`unmatched_subdicts:manufacturedBidomain` waiver. The remaining entries match
current scanner output.

If scanner quality is improved later, the useful change would be an explicit
path/access-mode model or a real C++ parser, together with a regression test
for duplicate leaf names such as myocardium `cm` versus Purkinje `cm`. That is
an audit-quality improvement, not a required solver correction.
