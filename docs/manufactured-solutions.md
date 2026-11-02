# Manufactured solutions

cardiacFoam verifies its monodomain and eikonal solvers against manufactured
solutions: exact fields, with the sources that make them exact, for which the
solver error can be measured directly. Each manufactured solution is defined
once, as a sympy generator, and everything else is produced from it.

## Where things are

| You want to | Go to |
|---|---|
| read how each manufactured solution is built | [monodomain_mms.pdf](../tutorials/manufacturedSolutions/generators/docs/monodomain/monodomain_mms.pdf), [eikonal_mms.pdf](../tutorials/manufacturedSolutions/generators/docs/eikonal/eikonal_mms.pdf) |
| use an exact solution in your own code | [`generators/reference/`](../tutorials/manufacturedSolutions/generators/reference): one NumPy-only module per case |
| change or add a manufactured solution | [`generators/`](../tutorials/manufacturedSolutions/generators/README.md): the sympy generators and their self-checks |
| create a new manufactured solution | "Creating a new manufactured solution" in [`generators/README.md`](../tutorials/manufacturedSolutions/generators/README.md) and in both PDFs, starting from `template_mms.py` |
| run a verification case in cardiacFoam | [`tutorials/manufacturedSolutions/`](../tutorials/manufacturedSolutions/README.md): one folder per case, with `Allrun` and a regression |
| see how the solver compares with the exact solution | [`src/verificationModels/`](../src/verificationModels/README.md): the verifiers and the C++ reference headers |

The [`generators/`](../tutorials/manufacturedSolutions/generators/README.md)
folder stands alone: it needs Python only (plus Quarto and TeX to rebuild the
PDFs), not cardiacFoam or OpenFOAM, so it can be copied and used by any code.

## Cases

| Model | Case | Generator |
|---|---|---|
| monodomain | FDA cube, d = 1, 2, 3 | `fda_monodomain_mms.py` |
| monodomain | rotated anisotropic cube | `rotated_monodomain_mms.py` |
| monodomain | sealed spherical shell | `shell_monodomain_mms.py` |
| monodomain | pseudo-ECG reference | `pseudo_ecg_mms.py` |
| monodomain | Purkinje graph | `graph_monodomain_mms.py` |
| monodomain | 1D–3D Purkinje–myocardium coupling | `coupled_1d3d_mms.py` |
| eikonal | Cartesian cube, d = 1, 2, 3 | `cartesian_eikonal_mms.py` |
| eikonal | eikonal-driven ECG reference | `eikonal_ecg_mms.py` |
| eikonal | sealed spherical shell, any `c0` | `shell_eikonal_mms.py` |

## How the pieces connect

```text
generators/*_mms.py  (sympy, checked before writing)
   ├── reference/*.py        NumPy exact solutions, for any code
   ├── latex/*.tex           formulas printed by the two PDFs
   └── src/verificationModels/**/*Reference.H   C++ used by the solver
```

The C++ headers are generated, whole or between
`// BEGIN GENERATED <name>` and `// END GENERATED <name>` markers, and
`python <generator> --check <header>` fails if a header no longer matches its
generator. `generators/parity/` compiles the headers and compares them with the
NumPy references. Edit the generator, not the header or the formulas in a
README.

## Quick start

```bash
cd tutorials/manufacturedSolutions/generators
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python shell_eikonal_mms.py --self-check --c0 5
.venv/bin/python -c "from reference import fda_monodomain as m; print(m.Vm_d3(0.2, 0.3, 0.4, 0.1))"
```

Bidomain, bath-bidomain and electromechanics manufactured solutions are not
yet in this scheme; their references are in `src/verificationModels/` and
their cases in `tutorials/manufacturedSolutions/`.
