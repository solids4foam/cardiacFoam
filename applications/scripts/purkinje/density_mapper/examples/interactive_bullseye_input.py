from purkinje_density.interactive_bullseye import (
    DEFAULT_AHA17_REFERENCE_CFG,
    collect_segment_percentages_interactive,
)


def main() -> None:
    values = collect_segment_percentages_interactive(
        reference_cfg=DEFAULT_AHA17_REFERENCE_CFG,
        title="Select Segment Percentages (AHA17)",
        prompt_missing_after_close=True,
    )
    print("Final segment percentage map:")
    for seg_id, value in values.items():
        print(f"  {seg_id}: {value:.1f}%")


if __name__ == "__main__":
    main()
