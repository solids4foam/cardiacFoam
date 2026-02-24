"""Pipeline package and internal step registration."""

from cardiac_preproc.pipeline.registry import register_step


def register_default_steps() -> None:
    from cardiac_preproc.steps.diffusivity import run_from_context as run_diffusivity_step
    from cardiac_preproc.steps.purkinje_slab import (
        run_from_context as run_purkinje_slab_step,
    )
    from cardiac_preproc.steps.scar import run_from_context as run_scar_step

    register_step("diffusivity", run_diffusivity_step)
    register_step("purkinje_slab", run_purkinje_slab_step)
    register_step("scar", run_scar_step)
