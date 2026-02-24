from purkinje_density.tree_metrics import (
    TreeSample,
    compute_tree_presence_metrics,
)


def main() -> None:
    samples = [
        TreeSample(tissue="endocardium", is_terminal=False),
        TreeSample(tissue="myocardium", is_terminal=False),
        TreeSample(tissue="endocardium", is_terminal=False),
        TreeSample(tissue="endocardium", is_terminal=True),
        TreeSample(tissue="myocardium", is_terminal=True),
        TreeSample(tissue="endocardium", is_terminal=True),
    ]

    metrics = compute_tree_presence_metrics(samples)
    print("P(any part in tissue):", metrics.any_part_probability_by_tissue)
    print("P(terminal end in tissue):", metrics.terminal_probability_by_tissue)


if __name__ == "__main__":
    main()
