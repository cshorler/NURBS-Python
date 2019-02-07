"""
    Tests for the NURBS-Python package
    Released under The MIT License. See LICENSE file for details.
    Copyright (c) 2018-2019 Onur Rauf Bingol

    Requires "pytest" to run.
"""

from pytest import fixture, mark
from geomdl import BSpline
from geomdl import evaluators
from geomdl import convert

GEOMDL_DELTA = 0.001


@fixture
def spline_surf():
    """ Creates a B-spline surface instance """
    # Create a surface instance
    surf = BSpline.Surface()

    # Set degrees
    surf.degree_u = 3
    surf.degree_v = 3

    ctrlpts = [[-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0],
               [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0], [-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0],
               [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0],
               [-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0],
               [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0], [5.0, -25.0, -3.0], [5.0, -15.0, -2.0],
               [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0],
               [15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0],
               [15.0, 15.0, -4.0], [15.0, 25.0, -8.0], [25.0, -25.0, -10.0], [25.0, -15.0, -5.0],
               [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]]

    # Set control points
    surf.set_ctrlpts(ctrlpts, 6, 6)

    # Set knot vectors
    surf.knotvector_u = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]
    surf.knotvector_v = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]

    return surf


def test_bspline_curve_name(spline_surf):
    spline_surf.name = "Surface Testing"
    assert spline_surf.name == "Surface Testing"


def test_bspline_surface_degree_u(spline_surf):
    assert spline_surf.degree_u == 3


def test_bspline_surface_degree_v(spline_surf):
    assert spline_surf.degree_v == 3


def test_bspline_surface_ctrlpts(spline_surf):
    assert spline_surf.ctrlpts2d[1][1] == [-15.0, -15.0, -4.0]
    assert spline_surf.dimension == 3


def test_bspline_surface_knot_vector_u(spline_surf):
    assert spline_surf.knotvector_u == [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]


def test_bspline_surface_knot_vector_v(spline_surf):
    assert spline_surf.knotvector_v == [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]


@mark.parametrize("param, res", [
    ((0.0, 0.0), (-25.0, -25.0, -10.0)),
    ((0.0, 0.2), (-25.0, -11.403, -3.385)),
    ((0.0, 1.0), (-25.0, 25.0, -10.0)),
    ((0.3, 0.0), (-7.006, -25.0, -5.725)),
    ((0.3, 0.4), [-7.006, -3.308, -6.265]),
    ((0.3, 1.0), [-7.006, 25.0, -5.725]),
    ((0.6, 0.0), (3.533, -25.0, -4.224)),
    ((0.6, 0.6), (3.533, 3.533, -6.801)),
    ((0.6, 1.0), (3.533, 25.0, -4.224)),
    ((1.0, 0.0), (25.0, -25.0, -10.0)),
    ((1.0, 0.8), (25.0, 11.636, -2.751)),
    ((1.0, 1.0), (25.0, 25.0, -10.0))
])
def test_bspline_surface_eval(spline_surf, param, res):
    evalpt = spline_surf.evaluate_single(param)
    assert abs(evalpt[0] - res[0]) < GEOMDL_DELTA
    assert abs(evalpt[1] - res[1]) < GEOMDL_DELTA
    assert abs(evalpt[2] - res[2]) < GEOMDL_DELTA


def test_bspline_surface_deriv(spline_surf):
    der1 = spline_surf.derivatives(u=0.35, v=0.35, order=2)
    spline_surf.evaluator = evaluators.SurfaceEvaluator2()
    der2 = spline_surf.derivatives(u=0.35, v=0.35, order=2)
    for k in range(0, 3):
        for l in range(0, 3 - k):
            assert abs(der1[k][l][0] - der2[k][l][0]) < GEOMDL_DELTA
            assert abs(der1[k][l][1] - der2[k][l][1]) < GEOMDL_DELTA
            assert abs(der1[k][l][2] - der2[k][l][2]) < GEOMDL_DELTA


@mark.parametrize("params, uv, res", [
    (dict(u=0.3, v=0.4), (0.3, 0.4), (-7.006, -3.308, -6.265)),
    (dict(u=0.3, num_u=2), (0.3, 0.4), (-7.006, -3.308, -6.265)),
    (dict(v=0.3, num_v=2), (0.3, 0.4), (-7.006, -3.308, -6.265))
])
def test_bspline_surface_insert_knot_eval(spline_surf, params, uv, res):
    # Insert knot
    spline_surf.insert_knot(**params)

    # Evaluate surface
    evalpt = spline_surf.evaluate_single(uv)

    assert abs(evalpt[0] - res[0]) < GEOMDL_DELTA
    assert abs(evalpt[1] - res[1]) < GEOMDL_DELTA
    assert abs(evalpt[2] - res[2]) < GEOMDL_DELTA


@mark.parametrize("params, idx, val", [
    (dict(v=0.3, num_v=2), 4, 0.3),
    (dict(v=0.3, num_v=2), 6, 0.33)
])
def test_bspline_surface_insert_knot_kv_v(spline_surf, params, idx, val):
    # Insert knot
    spline_surf.insert_knot(**params)

    assert spline_surf.knotvector_v[idx] == val


@mark.parametrize("params, idx, val", [
    (dict(u=0.33, num_u=2), 3, 0.0),
    (dict(u=0.33, num_u=1), 6, 0.66)
])
def test_bspline_surface_insert_kv_u(spline_surf, params, idx, val):
    # Insert knot
    spline_surf.insert_knot(**params)

    assert spline_surf.knotvector_u[idx] == val


@fixture
def nurbs_surf(spline_surf):
    surf = convert.bspline_to_nurbs(spline_surf)
    return surf


def test_nurbs_weights(nurbs_surf):
    assert len(nurbs_surf.weights) == nurbs_surf.ctrlpts_size
    assert nurbs_surf.weights[5] == 1.0


@mark.parametrize("param, res", [
    ((0.0, 0.0), (-25.0, -25.0, -10.0)),
    ((0.0, 0.2), (-25.0, -11.403, -3.385)),
    ((0.0, 1.0), (-25.0, 25.0, -10.0)),
    ((0.3, 0.0), (-7.006, -25.0, -5.725)),
    ((0.3, 0.4), [-7.006, -3.308, -6.265]),
    ((0.3, 1.0), [-7.006, 25.0, -5.725]),
    ((0.6, 0.0), (3.533, -25.0, -4.224)),
    ((0.6, 0.6), (3.533, 3.533, -6.801)),
    ((0.6, 1.0), (3.533, 25.0, -4.224)),
    ((1.0, 0.0), (25.0, -25.0, -10.0)),
    ((1.0, 0.8), (25.0, 11.636, -2.751)),
    ((1.0, 1.0), (25.0, 25.0, -10.0))
])
def test_nurbs_surface_eval(nurbs_surf, param, res):
    evalpt = nurbs_surf.evaluate_single(param)
    assert abs(evalpt[0] - res[0]) < GEOMDL_DELTA
    assert abs(evalpt[1] - res[1]) < GEOMDL_DELTA
    assert abs(evalpt[2] - res[2]) < GEOMDL_DELTA


@mark.parametrize("param, order, res", [
    ((0.0, 0.25), 1, [[[10.4284, 3.7864], [76.3058, 9.4804]], [[-103.5570, -23.8313], [0.0, 0.0, 0.0]]]),
    ((0.95, 0.75), 2, [[[-8.41, -3.7471], [50.8270, 5.0301], [-363.4711, 10.6923]], [[-91.7390, -27.204], [1460.1162, 440.0493], [0.0, 0.0, 0.0]], [[-1738.5193, -535.8413], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]])
])
def test_nurbs_surface_deriv(nurbs_surf, param, order, res):
    deriv = nurbs_surf.derivatives(*param, order=order)

    for computed, expected in zip(deriv, res):
        for idx in range(order + 1):
            for c, e in zip(computed[idx], expected[idx]):
                assert abs(c - e) < GEOMDL_DELTA


def test_surface_bounding_box(spline_surf):
    # Evaluate bounding box
    to_check = spline_surf.bbox

    # Evaluation result
    result = ((-25.0, -25.0, -10.0), (25.0, 25.0, 2.0))

    assert abs(to_check[0][0] - result[0][0]) < GEOMDL_DELTA
    assert abs(to_check[0][1] - result[0][1]) < GEOMDL_DELTA
    assert abs(to_check[0][2] - result[0][2]) < GEOMDL_DELTA
    assert abs(to_check[1][0] - result[1][0]) < GEOMDL_DELTA
    assert abs(to_check[1][1] - result[1][1]) < GEOMDL_DELTA
    assert abs(to_check[1][2] - result[1][2]) < GEOMDL_DELTA
