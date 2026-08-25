import inspect

from openfoam_driver.plugins.cardiacfoam.tutorials.registry import SPEC_FACTORIES


_DEAD_POSTPROCESS_SELECTOR_NAMES = {
    "cv_extract_script_relpath",
    "postprocess_function_name",
    "postprocess_script_relpath",
    "table_summary_relpath",
}


def test_registered_tutorials_do_not_advertise_undiscoverable_postprocess_selectors():
    """Postprocessing discovers setup scripts by their standard entry point."""
    checked_factories = set()
    for factory in SPEC_FACTORIES.values():
        if factory in checked_factories:
            continue
        checked_factories.add(factory)
        advertised = set(inspect.signature(factory).parameters)
        assert advertised.isdisjoint(_DEAD_POSTPROCESS_SELECTOR_NAMES), (
            f"{factory.__module__}.{factory.__name__} advertises a postprocess "
            "selector that the runtime cannot discover"
        )
