# frontend

This directory contains the Next.js frontend prototype for the local
`cardiacFoam` driver/UI workflow.

## Intended role

The frontend is meant to consume the machine-readable metadata exposed by:

- `foamctl describe`
- `run_manifest.json`
- `plots.json`

and provide a GUI around tutorial selection, parameter editing, and run/result
inspection.

## Current status

This is still a prototype workspace, not the authoritative architecture source.
The backend contract lives primarily in:

- `applications/scripts/driverFoam/openfoam_driver/README.md`
- `applications/scripts/driverFoam/openfoam_driver/GUI_CONTRACT.md`
- `applications/scripts/driverFoam/openfoam_driver/FRONTEND_HANDOFF.md`

## Development

```bash
cd frontend
npm install
npm run dev
```
