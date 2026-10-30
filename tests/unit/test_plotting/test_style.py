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


# --- Neighbouring panel labels stay tellable apart ---------------------------

# Panel labels whose clusters touch in real PC space: the European/Middle-
# Eastern block, the African pair, the Asian groups, and the admixed AMR and
# CAH against every group they spread towards.
NEIGHBOURS = [
    ("EUR", "FIN"), ("EUR", "AJ"), ("EUR", "MDE"), ("AJ", "MDE"), ("FIN", "AJ"),
    ("FIN", "MDE"), ("AFR", "AAC"), ("CAS", "EUR"), ("CAS", "SAS"), ("CAS", "EAS"),
    ("CAS", "MDE"), ("CAS", "AMR"), ("SAS", "EAS"), ("SAS", "MDE"), ("AMR", "EUR"),
    ("AMR", "EAS"), ("AMR", "AFR"), ("AMR", "AAC"), ("AMR", "SAS"), ("AMR", "MDE"),
    ("AMR", "FIN"), ("CAH", "AFR"), ("CAH", "EUR"), ("CAH", "AAC"), ("CAH", "AMR"),
    ("CAH", "MDE"), ("AJ", "AFR"), ("MDE", "AFR"), ("MDE", "AAC"),
]

# Machado et al. (2009) full-severity deuteranopia and protanopia, linear RGB.
_DEUTAN = [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
           [-0.011820, 0.042940, 0.968881]]
_PROTAN = [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216],
           [-0.003882, -0.048116, 1.051998]]


def _linear(hex_color: str):
    import numpy as np

    c = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)]) / 255
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _oklab(rgb):
    import numpy as np

    m1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                   [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    m2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                   [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    return m2 @ np.cbrt(m1 @ rgb)


def _delta_e(a: str, b: str, sim=None) -> float:
    """OKLab distance x100, optionally under a colour-vision simulation."""
    import numpy as np

    la, lb = _linear(a), _linear(b)
    if sim is not None:
        la, lb = (np.clip(np.array(sim) @ x, 0, 1) for x in (la, lb))
    return float(np.linalg.norm(_oklab(la) - _oklab(lb)) * 100)


def test_neighbouring_panel_labels_clear_both_floors() -> None:
    """The dataviz floors: normal-vision dE >= 15, CVD dE >= 8."""
    weak = []
    for a, b in NEIGHBOURS:
        ca, cb = PANEL_COLORS[a], PANEL_COLORS[b]
        normal = _delta_e(ca, cb)
        cvd = min(_delta_e(ca, cb, _DEUTAN), _delta_e(ca, cb, _PROTAN))
        if normal < 15 or cvd < 8:
            weak.append(f"{a}-{b}: normal {normal:.1f}, cvd {cvd:.1f}")
    assert not weak, weak
