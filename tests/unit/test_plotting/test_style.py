"""Shared colours and layout helpers."""

from genotools.plotting.style import (
    CATEGORICAL,
    MARKERS,
    NEUTRAL,
    PANEL_COLORS,
    PANEL_MARKERS,
    grid,
    label_colors,
    label_styles,
    step_label,
)


class TestLabelStyles:
    def test_panel_labels_keep_their_colour(self) -> None:
        """EUR is the same colour whichever other groups a run has."""
        assert label_colors(["EUR"])["EUR"] == PANEL_COLORS["EUR"]
        assert label_colors(list(PANEL_COLORS))["EUR"] == PANEL_COLORS["EUR"]

    def test_panel_colours_are_distinct_slots(self) -> None:
        assert len(set(PANEL_COLORS.values())) == len(PANEL_COLORS)
        assert set(PANEL_COLORS.values()) <= set(CATEGORICAL)

    def test_order_of_input_does_not_matter(self) -> None:
        assert label_styles(["SAS", "X1", "EUR"]) == label_styles(["EUR", "SAS", "X1"])

    def test_custom_vocabulary_gets_spare_slots(self) -> None:
        """A --ref-labels vocabulary the panel never saw still gets colours."""
        colors = label_colors(["Pop1", "Pop2", "EUR"])
        assert colors["EUR"] == PANEL_COLORS["EUR"]
        assert colors["Pop1"] != colors["Pop2"]
        assert {colors["Pop1"], colors["Pop2"]} <= set(CATEGORICAL) - {PANEL_COLORS["EUR"]}

    def test_no_generated_hue_past_the_slots(self) -> None:
        labels = [f"L{i:02d}" for i in range(14)]
        colors = label_colors(labels)
        assert [colors[l] for l in labels[: len(CATEGORICAL)]] == CATEGORICAL
        assert {colors[l] for l in labels[len(CATEGORICAL):]} == {NEUTRAL}

    def test_every_panel_label_has_its_own_shape(self) -> None:
        styles = label_styles(PANEL_COLORS)
        assert len({s.marker for s in styles.values()}) == len(PANEL_COLORS)
        assert len({s.symbol for s in styles.values()}) == len(PANEL_COLORS)
        assert all((s.marker, s.symbol) in MARKERS for s in styles.values())
        assert set(PANEL_MARKERS) == set(PANEL_COLORS)

    def test_custom_labels_avoid_the_panel_shapes_present(self) -> None:
        styles = label_styles(["EUR", "Pop1"])
        assert styles["EUR"].marker == PANEL_MARKERS["EUR"][0]
        assert styles["Pop1"].marker != styles["EUR"].marker

    def test_labels_are_strings(self) -> None:
        assert set(label_colors([1, "1"])) == {"1"}


def test_step_label() -> None:
    assert step_label("het_prune") == "het"
    assert step_label("case_control_missingness_prune") == "case_control"


def test_grid() -> None:
    assert grid(1) == (1, 1)
    assert grid(4) == (2, 3)
    assert grid(0) == (0, 1)
