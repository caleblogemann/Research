from sageresearch.discontinuousgalerkin import canonical_element
from sageresearch.utils import symbolic_vector_matrix as svm

import sage.all as sa

one = sa.Integer(1)
two = sa.Integer(2)

def get_q_variables(num_eqns):
    return svm.get_vector_variable('q', num_eqns)


def get_t_variable():
    return sa.SR.symbol('t')


def get_advection_flux_function(wavespeed):
    q = get_q_variables(1)
    x = canonical_element.get_mesh_variables_1d()
    t = get_t_variable()
    f = (wavespeed * q[0]).function(q[0], x, t)
    return f


def get_burgers_flux_function():
    q = get_q_variables(1)
    x = canonical_element.get_mesh_variables_1d()
    t = get_t_variable()
    f = (one / two * q[0] ** two).function(q[0], x, t)
    return f
