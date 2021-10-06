import sage.all as sa


def _get_legendre_polynomials_on_interval(
    max_order, variable, lower_bound=-1.0, upper_bound=1.0
):
    # polynomials orthogonal over (lower_bound, upper_bound)
    # not normalized in any way

    # linear_transformation from (-1, 1) to (lower_bound, upper_bound)
    delta_x = upper_bound - lower_bound
    assert delta_x > 0.0
    x_c = (upper_bound + lower_bound) / sa.Integer(2)
    dict_ = {variable: (variable - x_c) * sa.Integer(2) / delta_x}
    # legendre_P give orthogonal polynomials over (-1, 1)
    legendre_array = [
        sa.legendre_P(i, variable).function(variable).subs(dict_)
        for i in range(max_order)
    ]

    return legendre_array


def get_legendre_polynomials(
    max_order,
    normalization_constant,
    variable,
    lower_bound=-sa.Integer(1),
    upper_bound=sa.Integer(1),
):
    # polynomials orthonormal over (lower_bound, upper_bound)
    # with normalization_constant given
    # (normalization_constant * phi[i] * phi[j]).integrate(x, lower_bound, upper_bound)
    #   = \delta_{ij}

    legendre_array = _get_legendre_polynomials_on_interval(
        max_order, variable, lower_bound, upper_bound
    )

    # need to normalize
    norms_squared = [
        (legendre_array[l] * legendre_array[l]).integrate(
            variable, -sa.Integer(1), sa.Integer(1)
        )
        for l in range(max_order)
    ]
    phi = [
        (
            sa.Integer(1)
            / sa.sqrt(norms_squared[l] * normalization_constant)
            * legendre_array[l]
        )
        for l in range(max_order)
    ]
    return phi


def get_legendre_polynomials_fixed_lower_endpoint(
    max_order,
    variable,
    fixed_endpoint=sa.Integer(1),
    lower_bound=-sa.Integer(1),
    upper_bound=sa.Integer(1),
):
    # polynomials orthogonal over an interval (lower_bound, upper_bound)
    # such that phi(lower_bound) = fixed_endpoint

    legendre_array = _get_legendre_polynomials_on_interval(
        max_order, variable, lower_bound, upper_bound
    )

    # scale so phi(lower_bound) = fixed_endpoint
    phi = [
        (
            legendre_array[i] * fixed_endpoint / legendre_array[i](lower_bound)
        ).full_simplify()
        for i in range(max_order)
    ]
    return phi


def check_orthogonality(
    phi, variable, lower_bound=-sa.Integer(1), upper_bound=sa.Integer(1)
):
    max_order = len(phi)
    for l in range(max_order):
        for m in range(max_order):
            actual_integral = (phi[l] * phi[m]).integrate(
                variable, lower_bound, upper_bound
            )
            if l != m:
                assert actual_integral <= 1e-12
    return True


def check_normality(
    phi,
    normalization_constant,
    variable,
    lower_bound=-sa.Integer(1),
    upper_bound=sa.Integer(1),
):
    max_order = len(phi)
    for l in range(max_order):
        actual_integral = (normalization_constant * phi[l] * phi[l]).integrate(
            variable, lower_bound, upper_bound
        )
        assert abs(actual_integral - 1.0) <= 1e-12
    return True


def get_mass_matrix(
    phi, variable, lower_bound=-sa.Integer(-1), upper_bound=sa.Integer(1)
):
    max_order = len(phi)
    mass_matrix = sa.matrix.identity(sa.QQ, max_order)
    for i in range(max_order):
        mass_matrix[i, i] = (phi[i] * phi[i]).integrate(
            variable, lower_bound, upper_bound
        )
    return mass_matrix


def get_mass_matrix_inverse(
    phi, variable, lower_bound=-sa.Integer(1), upper_bound=sa.Integer(1)
):
    max_order = len(phi)
    mass_matrix_inverse = sa.matrix.identity(sa.QQ, max_order)
    for i in range(max_order):
        mass_matrix_inverse[i, i] = sa.Integer(1) / (phi[i] * phi[i]).integrate(
            variable, lower_bound, upper_bound
        )
    return mass_matrix_inverse


if __name__ == "__main__":
    pass
