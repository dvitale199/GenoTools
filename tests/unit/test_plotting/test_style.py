"""Shared colours and layout helpers."""

from genotools.plotting.style import CATEGORICAL, NEUTRAL, grid, label_colors, step_label


class TestLabelColors:
    def test_colour_follows_the_label_not_its_rank(self) -> None:
        """A plot showing a subset must not repaint the survivors."""
        full = label_colors(["SAS", "EUR", "AFR", "AMR"])
        assert full == label_colors(["AMR", "AFR", "SAS", "EUR"])
        assert full["AFR"] == CATEGORICAL[0]
        assert full["SAS"] == CATEGORICAL[3]

    def test_no_ninth_hue(self) -> None:
        labels = [f"L{i:02d}" for i in range(11)]
        colors = label_colors(labels)
        assert [colors[l] for l in labels[:8]] == CATEGORICAL
        assert {colors[l] for l in labels[8:]} == {NEUTRAL}

    def test_labels_are_strings(self) -> None:
        assert set(label_colors([1, "1"])) == {"1"}


def test_step_label() -> None:
    assert step_label("het_prune") == "het"
    assert step_label("case_control_missingness_prune") == "case_control"


def test_grid() -> None:
    assert grid(1) == (1, 1)
    assert grid(4) == (2, 3)
    assert grid(0) == (0, 1)
