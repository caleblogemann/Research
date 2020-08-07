def _get_legendre_polynomials_on_interval(max_order, lower_bound=-1, upper_bound=1):
    # polynomials orthogonal over (lower_bound, upper_bound) 
    # not normalized in any way

    # linear_transformation from (-1, 1) to (lower_bound, upper_bound)
    delta_x = upper_bound - lower_bound
    assert delta_x > 0
    x_c = (upper_bound + lower_bound) / 2
    dict_ = {x: (x - x_c) * 2 / delta_x}
    # legendre_P give orthogonal polynomials over (-1, 1)
    legendre_array = [legendre_P(i, x).function(x).subs(dict_) for i in range(max_order)]

    return legendre_array


def get_legendre_polynomials(max_order, normalization_constant, lower_bound=-1, upper_bound=1):
    # polynomials orthonormal over (lower_bound, upper_bound) 
    # with normalization_constant given
    # (normalization_constant * phi[i] * phi[j]).integrate(x, lower_bound, upper_bound) = \delta_{ij}

    legendre_array = _get_legendre_polynomials_on_interval(max_order, lower_bound, upper_bound)

    # need to normalize
    norms_squared = [
        (legendre_array[l] * legendre_array[l]).integrate(x, -1, 1)
        for l in range(max_order)
    ]
    phi = [
        (1 / sqrt(norms_squared[l] * normalization_constant) * legendre_array[l])
        for l in range(max_order)
    ]
    return phi


def get_legendre_polynomials_fixed_lower_endpoint(max_order, fixed_endpoint=1, lower_bound=-1, upper_bound=1):
    # polynomials orthogonal over an interval (lower_bound, upper_bound)
    # such that phi(lower_bound) = fixed_endpoint

    legendre_array = _get_legendre_polynomials_on_interval(max_order, lower_bound, upper_bound)

    # scale so phi(lower_bound) = fixed_endpoint
    phi = [(legendre_array[i] * fixed_endpoint / legendre_array[i](lower_bound)).full_simplify() for i in range(max_order)]
    return phi


def check_orthogonality(phi, lower_bound=-1, upper_bound=1):
    max_order = len(phi)
    for l in range(max_order):
        for m in range(max_order):
            actual_integral = (phi[l] * phi[m]).integrate(x, lower_bound, upper_bound)
            if l != m:
                assert actual_integral == 0
    return True


def check_normality(phi, normalization_constant, lower_bound=-1, upper_bound=1):
    max_order = len(phi)
    for l in range(max_order):
        actual_integral = (normalization_constant * phi[l] * phi[l]).integrate(x, lower_bound, upper_bound)
        assert actual_integral == 1
    return True


def get_mass_matrix(phi, lower_bound=-1, upper_bound=1):
    max_order = len(phi)
    mass_matrix = matrix.identity(QQ, max_order)
    for i in range(max_order):
        mass_matrix[i, i] = (phi[i] * phi[i]).integrate(x, lower_bound, upper_bound)
    return mass_matrix


def get_mass_matrix_inverse(max_order, lower_bound=-1, upper_bound=1):
    max_order = len(phi)
    mass_matrix_inverse = matrix.identity(QQ, max_order)
    for i in range(max_order):
        mass_matrix_inverse[i, i] = 1 / (phi[i] * phi[i]).integrate(x, lower_bound, upper_bound)
    return mass_matrix_inverse


if __name__ == "__main__":
    max_order = 4
    normalization_constant = 1 / 2
    phi = get_legendre_polynomials(max_order, normalization_constant)
    check_orthogonality(phi)
    check_normality(phi, normalization_constant)
    m = get_mass_matrix(phi)
    m_inv = get_mass_matrix_inverse(phi)
    assert (m * m_inv - matrix.identity(max_order)).norm() == 0

    max_order = 10
    fixed_endpoint = 1
    lower_bound = 0
    upper_bound = 1
    phi = get_legendre_polynomials_fixed_lower_endpoint(max_order, fixed_endpoint, lower_bound, upper_bound)
    check_orthogonality(phi, lower_bound, upper_bound)
    m = get_mass_matrix(phi, lower_bound, upper_bound)
    m_inv = get_mass_matrix_inverse(phi, lower_bound, upper_bound)
    assert (m * m_inv - matrix.identity(max_order)).norm() == 0
    print("LegendrePolynomials.sage tests passed")
