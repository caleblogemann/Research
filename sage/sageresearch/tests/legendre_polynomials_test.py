import sageresearch.utils.legendre_polynomials as legendre_polynomials

import sage.all as sa


def test_get_legendre_polynomials():
    max_order = 4
    normalization_constant = 0.5
    phi = legendre_polynomials.get_legendre_polynomials(
        max_order, normalization_constant
    )

    legendre_polynomials.check_orthogonality(phi)
    legendre_polynomials.check_normality(phi, normalization_constant)
    m = legendre_polynomials.get_mass_matrix(phi)
    m_inv = legendre_polynomials.get_mass_matrix_inverse(phi)
    assert (m * m_inv - sa.matrix.identity(max_order)).norm() <= 1e-12


def test_get_legendre_polynomials_fixed_lower_endpoint():
    max_order = 10
    fixed_endpoint = 1.0
    lower_bound = 0.0
    upper_bound = 1.0
    phi = legendre_polynomials.get_legendre_polynomials_fixed_lower_endpoint(
        max_order, fixed_endpoint, lower_bound, upper_bound
    )
    legendre_polynomials.check_orthogonality(phi, lower_bound, upper_bound)
    m = legendre_polynomials.get_mass_matrix(phi, lower_bound, upper_bound)
    m_inv = legendre_polynomials.get_mass_matrix_inverse(phi, lower_bound, upper_bound)
    assert (m * m_inv - sa.matrix.identity(max_order)).norm() <= 1e-12
