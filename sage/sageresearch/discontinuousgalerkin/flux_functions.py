from sageresearch.discontinuousgalerkin import canonical_element
from sageresearch.utils import symbolic_vector_matrix as svm

import sage.all as sa

one = sa.Integer(1)
two = sa.Integer(2)


def get_q_variables(num_eqns):
    return svm.get_vector_variable("q", num_eqns)


def get_t_variable():
    return sa.SR.symbol("t")


def get_advection_1d(wavespeed):
    q = get_q_variables(1)
    x = canonical_element.Interval.get_mesh_variables()
    t = get_t_variable()
    f = (wavespeed * q[0]).function(q[0], x, t)
    return f


def get_advection_2d(wavespeed):
    q = get_q_variables(1)
    tuple_ = canonical_element.CanonicalElement2D.get_mesh_variables()
    x = tuple_[0]
    y = tuple_[1]
    t = get_t_variable()
    f = sa.vector([wavespeed[0] * q[0], wavespeed[1] * q[0]]).function(q[0], x, y, t)
    return f


def get_burgers_1d():
    q = get_q_variables(1)
    x = canonical_element.Interval.get_mesh_variables()
    t = get_t_variable()
    f = (one / two * q[0] ** two).function(q[0], x, t)
    return f


def get_burgers_2d(n):
    # n should be unit vector
    q = get_q_variables()
    tuple_ = canonical_element.CanonicalElement2D.get_mesh_variables()
    x = tuple_[0]
    y = tuple_[1]
    t = get_t_variable()
    f = sa.vector(
        [n[0] * one / two * q[0] ** two, n[1] * one / two * q[0] ** two]
    ).function(q[0], x, y, t)
    return f


def get_shallow_water_1d():
    q = get_q_variables(3)
    pass


def get_shallow_water_2d():
    pass
