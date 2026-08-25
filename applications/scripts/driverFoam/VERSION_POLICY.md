# Version Policy

driverFOAM strictly adheres to [Semantic Versioning 2.0.0](https://semver.org/).

Given the nature of the project as a workflow orchestrator, we guarantee backward compatibility on the following components:

- **Command-line Interface (CLI)**: CLI arguments and standard exit codes.
- **Schemas**: The structural schema of `run-document.json`.
- **Plugin API**: The interface expected by `driverfoam.plugins` entry points.

Any breaking changes to the above will result in a major version bump. Internal Python modules not explicitly exposed as part of the public API may change in minor versions.
