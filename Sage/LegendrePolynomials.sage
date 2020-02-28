def get_legendre_polynomials(max_order, normalization_constant):
    # legendre_P give orthogonal polynomials over (-1, 1)
    legendre_array = [legendre_P(i, x).function(x) for i in range(max_order)]
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


def check_orthogonality(phi):
    max_order = len(phi)
    for l in range(max_order):
        for m in range(max_order):
            actual_integral = (phi[l] * phi[m]).integrate(x, -1, 1)
            if l != m:
                assert actual_integral == 0
    return True


def check_normality(phi, normalization_constant):
    max_order = len(phi)
    for l in range(max_order):
        actual_integral = (normalization_constant * phi[l] * phi[l]).integrate(x, -1, 1)
        assert actual_integral == 1
    return True


def get_mass_matrix(max_order, normalization_constant):
    return 1 / normalization_constant * matrix.identity(max_order)


def get_mass_matrix_inverse(max_order, normalization_constant):
    return normalization_constant * matrix.identity(max_order)


if __name__ == "__main__":
    phi = get_legendre_polynomials(10, 1 / 2)
    check_orthogonality(phi)
    check_normality(phi, 1 / 2)
    print("LegendrePolynomials.sage tests passed")
