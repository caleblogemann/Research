import sageresearch.shallowwater.hyperbolic_shallow_water_moment_equations as hswme
import sageresearch.shallowwater.shallow_water_moment_equations as swme

import sage.all as sa

max_num_moments = 3

def test_get_hswme_equations():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        A = hswme.get_hswme_equations(num_moments)
        assert A.nrows() == num_eqns
        assert A.ncols() == num_eqns

        # A should be equal to shallow water moment equations quasilinear matrix
        # with alpha_2 and higher set to zero
        tuple_ = swme.get_generalized_shallow_water_equations_1d(num_moments)
        A_s = tuple_[3]
        p = swme.get_primitive_variables_1d(num_moments)
        dict_ = {p[i]: 0 for i in range(3, num_eqns)}
        assert A_s.subs(dict_) == A


def test_get_beta_hswme_equations():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        A = hswme.get_beta_hswme_equations(num_moments)
        assert A.nrows() == num_eqns
        assert A.ncols() == num_eqns


def test_get_hswme_eigenvalues():
    assert False


def test_get_beta_hswme_eigenvalues():
    assert False
