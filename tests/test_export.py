import csv
import math
from collections import Counter

import ezdxf
import numpy as np
import pytest

from limerocket import export
from limerocket.cli import main
from limerocket.units import P_STANDARD_PA
from tests.conftest import EXAMPLES


def test_stl_wall_is_watertight_and_has_the_revolved_volume(lime1):
    chamber = lime1.chamber
    tris = export.wall_solid(chamber)
    edges = Counter()
    for tri in np.round(tris, 9):
        pts = [tuple(p) for p in tri]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edges[(pts[a], pts[b])] += 1
    assert all(edges[(b, a)] == n for (a, b), n in edges.items())

    signed = np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6
    xo, ro = export.outer_wall(chamber)
    x = np.linspace(chamber.x_m[0], chamber.x_m[-1], 20_000)
    ring = np.interp(x, xo, ro) ** 2 - np.interp(x, chamber.x_m, chamber.r_m) ** 2
    pappus = math.pi * np.trapezoid(ring, x) * 1e9
    polygon = export.STL_SEGMENTS / (2 * math.pi) * math.sin(2 * math.pi / export.STL_SEGMENTS)
    assert signed == pytest.approx(pappus * polygon, rel=0.02)


def test_dxf_holds_the_closed_wall_profile_in_millimetres(lime1, tmp_path):
    path = tmp_path / "c.dxf"
    export.write_dxf(path, lime1.chamber)
    doc = ezdxf.readfile(path)
    walls = [e for e in doc.modelspace() if e.dxf.layer == "WALL"]
    assert len(walls) == 1
    wall = walls[0]
    assert wall.is_closed
    points = np.array([(v.dxf.location.x, v.dxf.location.y) for v in wall.vertices])
    xo, _ = export.outer_wall(lime1.chamber)
    assert len(points) == lime1.chamber.x_m.size + xo.size
    assert points[:, 1].min() == pytest.approx(lime1.chamber.throat_radius_m * 1e3, abs=1e-4)


def test_eng_file_follows_the_rasp_format(lime1, tmp_path):
    path = tmp_path / "m.eng"
    export.write_eng(path, lime1)
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith(";")]
    header = lines[0].split()
    assert len(header) == 7
    assert float(header[1]) == pytest.approx(152)
    assert float(header[4]) == pytest.approx(lime1.burn.expelled_mass_kg, rel=1e-3)
    data = np.array([[float(v) for v in ln.split()] for ln in lines[1:]])
    assert np.all(np.diff(data[:, 0]) > 0)
    assert data[0, 0] > 0 and data[-1, 1] == 0.0
    impulse = np.trapezoid(data[:, 1], data[:, 0])
    assert impulse == pytest.approx(lime1.burn.total_impulse_ns(P_STANDARD_PA), rel=0.01)


def test_rocketpy_csv_is_two_numeric_columns(lime1, tmp_path):
    path = tmp_path / "t.csv"
    export.write_rocketpy_csv(path, lime1)
    rows = np.array([[float(v) for v in row] for row in csv.reader(path.open())])
    assert rows.shape[1] == 2
    assert np.trapezoid(rows[:, 1], rows[:, 0]) == pytest.approx(
        lime1.burn.total_impulse_ns(P_STANDARD_PA), rel=1e-3
    )


def test_cli_writes_the_report_and_every_export(tmp_path, capsys):
    assert main(["design", str(EXAMPLES / "lime-1.yaml"), "-o", str(tmp_path)]) == 0
    names = {p.name for p in tmp_path.iterdir()}
    assert {"report.html", "overview.svg", "lime-1.eng", "lime-1-chamber.stl"} <= names
    report = (tmp_path / "report.html").read_text()
    assert report.count("<svg") == 4
    assert "NASA CEA" in report
    assert "Apogee" in capsys.readouterr().out


def test_cli_reports_spec_errors_without_a_traceback(tmp_path, capsys):
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: x\nengine: {thrust_n: -5}\n")
    assert main(["design", str(bad), "-o", str(tmp_path / "out")]) == 2
    assert "engine.thrust_n: must be > 0" in capsys.readouterr().err
