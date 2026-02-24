from purkinje_density.tree_metrics import TreeSample, compute_tree_density_metrics


def main() -> None:
    samples = [
        # Non-terminal tree parts
        TreeSample(tissue="endocardium", is_terminal=False, section_tag="lv_apical_septal"),
        TreeSample(tissue="myocardium", is_terminal=False, section_tag="lv_mid_lateral"),
        TreeSample(tissue="endocardium", is_terminal=False, section_tag="lv_basal_anterior"),
        # Terminal ends
        TreeSample(tissue="endocardium", is_terminal=True, section_tag="lv_apical_septal"),
        TreeSample(tissue="endocardium", is_terminal=True, section_tag="lv_apical_septal"),
        TreeSample(tissue="myocardium", is_terminal=True, section_tag="lv_mid_lateral"),
    ]

    tissue_surface_areas_mm2 = {
        "endocardium": 120.0,
        "myocardium": 300.0,
    }
    section_surface_areas_mm2 = {
        "lv_apical_septal": 40.0,
        "lv_mid_lateral": 60.0,
    }

    metrics = compute_tree_density_metrics(
        samples=samples,
        tissue_surface_areas_mm2=tissue_surface_areas_mm2,
        section_surface_areas_mm2=section_surface_areas_mm2,
    )

    print("1) P(any part in tissue):", metrics.any_part_probability_by_tissue)
    print("2) P(end in tissue):", metrics.terminal_probability_by_tissue)
    print("3) End density by tissue [ends/mm^2]:", metrics.terminal_density_by_tissue_mm2)
    print("4) End density by section [ends/mm^2]:", metrics.terminal_density_by_section_mm2)


if __name__ == "__main__":
    main()
