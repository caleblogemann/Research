from sage.rings.number_field import number_field
from sageresearch.utils import legendre_polynomials
from sageresearch.discontinuousgalerkin import canonical_element

import numpy.polynomial.legendre as legendre
import sage.all as sa

zero = sa.Integer(0)
one = sa.Integer(1)
two = sa.Integer(2)
three = sa.Integer(3)
four = sa.Integer(4)
five = sa.Integer(5)
six = sa.Integer(6)
seven = sa.Integer(7)
eight = sa.Integer(8)
nine = sa.Integer(9)
eleven = sa.Integer(11)
fifteen = sa.Integer(15)
sixteen = sa.Integer(16)
eighteen = sa.Integer(18)
nineteen = sa.Integer(19)
twenty_one = sa.Integer(21)
twenty_seven = sa.Integer(27)
thirty_three = sa.Integer(33)
thirty_five = sa.Integer(35)
thirty_nine = sa.Integer(39)
forty_five = sa.Integer(45)
one_hundred_five = sa.Integer(105)
one_thirty_five = sa.Integer(135)


def get_nodal_basis_1d(nodes):
    points = [[x, 0] for x in nodes]
    xi = canonical_element.get_canonical_variables_1d()
    R = sa.PolynomialRing(sa.RR, xi)
    phi = []
    for i in range(len(nodes)):
        points[i][1] = 1
        phi.append((sa.SR(R.lagrange_polynomial(points))).function(xi))
        points[i][1] = 0
    return phi


def get_gauss_legendre_nodal_basis_1d(space_order):
    num_nodes = space_order
    phi = legendre.Legendre.basis(num_nodes)
    nodes = [sa.RR(x) for x in phi.roots()]
    return get_nodal_basis_1d(nodes)


def get_gauss_lobatto_nodal_basis_1d(space_order):
    num_nodes = space_order
    if num_nodes == 1:
        nodes = [zero]
    else:
        phi = legendre.Legendre.basis(num_nodes - 2)
        nodes = [sa.RR(x) for x in phi.roots()]
        nodes.append(one)
        nodes.insert(0, -one)

    return get_nodal_basis_1d(nodes)


def get_legendre_basis_1d(space_order, inner_product_constant=one / two):
    xi = canonical_element.get_canonical_variables_1d()
    return legendre_polynomials.get_legendre_polynomials(
        space_order, inner_product_constant, xi
    )


def get_mass_matrix(phi, integration_function):
    num_basis_cpts = len(phi)
    mass_matrix = sa.matrix(sa.SR, num_basis_cpts, num_basis_cpts)
    for i_basis_cpt in range(num_basis_cpts):
        for j_basis_cpt in range(i_basis_cpt + 1):
            f = phi[i_basis_cpt] * phi[j_basis_cpt]
            integral = integration_function(f)
            mass_matrix[i_basis_cpt, j_basis_cpt] = integral
            mass_matrix[j_basis_cpt, i_basis_cpt] = integral
    return mass_matrix


def get_mass_matrix_1d(phi):
    return get_mass_matrix(phi, canonical_element.integrate_over_canonical_element_1d)


def get_legendre_basis_2d_rectangle(space_order, inner_product_constant=one / four):
    tuple_ = canonical_element.get_canonical_variables_2d()
    xi = tuple_[0]
    eta = tuple_[1]
    phi_1d = legendre_polynomials.get_legendre_polynomials(
        space_order, sa.sqrt(inner_product_constant), xi
    )

    phi = []
    for i in range(space_order):
        for j in range(space_order - 1):
            phi.append((phi_1d[i](xi) * phi_1d[j](eta)).function(xi, eta))

    return phi


def get_mass_matrix_2d_rectangle(phi):
    return get_mass_matrix(
        phi, canonical_element.integrate_over_canonical_element_2d_rectangle
    )


def get_modal_basis_2d_triangle(space_order, inner_product_constant=one / two):
    coeffs = get_modal_basis_coefficients_2d_triangle(
        space_order, inner_product_constant
    )
    num_basis_cpts = coeffs.nrows()

    tuple_ = canonical_element.get_canonical_variables_2d()
    xi = tuple_[0]
    eta = tuple_[1]

    phi = []
    for i_basis_cpt in range(num_basis_cpts):
        result = 0
        i = 0
        for degree in range(space_order):
            for eta_degree in range(degree + 1):
                xi_degree = degree - eta_degree
                result += (
                    coeffs[i_basis_cpt, i] * (xi ** xi_degree) * (eta ** eta_degree)
                )
                i += 1
        phi.append(result.function(xi, eta))

    return phi


def get_modal_basis_coefficients_2d_triangle(
    space_order, inner_product_constant=one / two
):
    # coefficients of basis functions in terms on monomials
    # ordering 1, xi, eta, xi^2, xi * eta, eta^2, xi^3, xi^2 * eta, eta^2 * xi, ...
    # Done up to fourth degree monomials works for up to 5th order in space
    # coefficients are for polynomials orthonormal with inner product constant of
    # 0.5

    sqrt2 = sa.sqrt(two)
    sqrt3 = sa.sqrt(three)
    sqrt5 = sa.sqrt(five)
    sqrt6 = sa.sqrt(six)
    sqrt7 = sa.sqrt(seven)
    sqrt15 = sa.sqrt(fifteen)
    sqrt35 = sa.sqrt(thirty_five)

    coeffs = sa.matrix(sa.SR, 15, 15)
    coeffs[0, 0] = one

    coeffs[1, 0] = sqrt2 / two
    coeffs[1, 1] = sqrt2 * three / two

    coeffs[2, 0] = sqrt6 / two
    coeffs[2, 1] = sqrt6 / two
    coeffs[2, 2] = sqrt6

    coeffs[3, 0] = -sqrt3 / two
    coeffs[3, 1] = sqrt3
    coeffs[3, 3] = sqrt3 * five / two

    coeffs[4, 0] = nine / four
    coeffs[4, 1] = six
    coeffs[4, 2] = nine / two
    coeffs[4, 3] = fifteen / four
    coeffs[4, 4] = fifteen / two

    coeffs[5, 0] = sqrt15 / four
    coeffs[5, 1] = sqrt15
    coeffs[5, 2] = sqrt15 * three / two
    coeffs[5, 3] = sqrt15 / four
    coeffs[5, 4] = sqrt15 * three / two
    coeffs[5, 5] = sqrt15 * three / two

    coeffs[6, 0] = -three / four
    coeffs[6, 1] = -fifteen / four
    coeffs[6, 3] = fifteen / four
    coeffs[6, 6] = thirty_five / four

    coeffs[7, 0] = sqrt3 / four
    coeffs[7, 1] = sqrt3 * nineteen / four
    coeffs[7, 2] = sqrt3 / two
    coeffs[7, 3] = sqrt3 * thirty_nine / four
    coeffs[7, 4] = sqrt3 * nine
    coeffs[7, 6] = sqrt3 * twenty_one / four
    coeffs[7, 7] = sqrt3 * twenty_one / two

    coeffs[8, 0] = sqrt5 * five / four
    coeffs[8, 1] = sqrt5 * twenty_seven / four
    coeffs[8, 2] = sqrt5 * fifteen / two
    coeffs[8, 3] = sqrt5 * thirty_three / four
    coeffs[8, 4] = sqrt5 * eighteen
    coeffs[8, 5] = sqrt5 * fifteen / two
    coeffs[8, 6] = sqrt5 * seven / four
    coeffs[8, 7] = sqrt5 * twenty_one / two
    coeffs[8, 8] = sqrt5 * twenty_one / two

    coeffs[9, 0] = sqrt7 / four
    coeffs[9, 1] = sqrt7 * nine / four
    coeffs[9, 2] = sqrt7 * three
    coeffs[9, 3] = sqrt7 * nine / four
    coeffs[9, 4] = sqrt7 * nine
    coeffs[9, 5] = sqrt7 * fifteen / two
    coeffs[9, 6] = sqrt7 / four
    coeffs[9, 7] = sqrt7 * three
    coeffs[9, 8] = sqrt7 * fifteen / two
    coeffs[9, 9] = sqrt7 * five

    coeffs[10, 0] = sqrt5 * three / eight
    coeffs[10, 1] = sqrt5 * -three / two
    coeffs[10, 3] = sqrt5 * -twenty_one / four
    coeffs[10, 6] = sqrt5 * seven / two
    coeffs[10, 10] = sqrt5 * sa.Integer(63) / eight

    coeffs[11, 0] = sqrt15 / -two
    coeffs[11, 1] = sqrt15 / -two
    coeffs[11, 2] = sqrt15 * -one
    coeffs[11, 3] = sqrt15 * twenty_one / four
    coeffs[11, 6] = sqrt15 * twenty_one / two
    coeffs[11, 7] = sqrt15 * twenty_one / two
    coeffs[11, 10] = sqrt15 * twenty_one / four
    coeffs[11, 11] = sqrt15 * twenty_one / two

    coeffs[12, 0] = five / two
    coeffs[12, 1] = forty_five / two
    coeffs[12, 2] = fifteen
    coeffs[12, 3] = sa.Integer(255) / four
    coeffs[12, 4] = sa.Integer(90)
    coeffs[12, 5] = fifteen
    coeffs[12, 6] = sa.Integer(115) / two
    coeffs[12, 7] = sa.Integer(285) / two
    coeffs[12, 8] = sa.Integer(75)
    coeffs[12, 10] = forty_five / four
    coeffs[12, 11] = one_thirty_five / two
    coeffs[12, 12] = one_thirty_five / two

    coeffs[13, 0] = sqrt35 * seven / sixteen
    coeffs[13, 1] = sqrt35 * nine / two
    coeffs[13, 2] = sqrt35 * twenty_one / four
    coeffs[13, 3] = sqrt35 * nine
    coeffs[13, 4] = sqrt35 * forty_five / two
    coeffs[13, 5] = sqrt35 * one_hundred_five / eight
    coeffs[13, 6] = sqrt35 * eleven / two
    coeffs[13, 7] = sqrt35 * sa.Integer(51) / two
    coeffs[13, 8] = sqrt35 * sa.Integer(30)
    coeffs[13, 9] = sqrt35 * thirty_five / four
    coeffs[13, 10] = sqrt35 * nine / sixteen
    coeffs[13, 11] = sqrt35 * twenty_seven / four
    coeffs[13, 12] = sqrt35 * one_thirty_five / eight
    coeffs[13, 13] = sqrt35 * forty_five / four

    coeffs[14, 0] = sqrt5 * three / sixteen
    coeffs[14, 1] = sqrt5 * three
    coeffs[14, 2] = sqrt5 * fifteen / four
    coeffs[14, 3] = sqrt5 * twenty_seven / four
    coeffs[14, 4] = sqrt5 * forty_five / two
    coeffs[14, 5] = sqrt5 * one_thirty_five / eight
    coeffs[14, 6] = sqrt5 * three
    coeffs[14, 7] = sqrt5 * forty_five / two
    coeffs[14, 8] = sqrt5 * forty_five
    coeffs[14, 9] = sqrt5 * one_hundred_five / four
    coeffs[14, 10] = sqrt5 * three / sixteen
    coeffs[14, 11] = sqrt5 * fifteen / four
    coeffs[14, 12] = sqrt5 * one_thirty_five / eight
    coeffs[14, 13] = sqrt5 * one_hundred_five / four
    coeffs[14, 14] = sqrt5 * one_hundred_five / eight

    coeffs = coeffs / sa.sqrt(canonical_element.volume_2d_triangle * inner_product_constant)
    num_basis_cpts = int(space_order * (space_order + 1) / 2)
    return coeffs[:num_basis_cpts, :num_basis_cpts]


def compute_modal_basis_coefficients_2d_triangle(space_order):
    tuple_ = canonical_element.get_canonical_variables_2d()
    xi = tuple_[0]
    eta = tuple_[1]
    R.<xi, eta> = sa.PolynomialRing(sa.SR, order='negdeglex')

    def integral_over_triangle(poly):
        print(type(poly))
        xi_integral = poly.integral(xi)
        print(type(xi_integral))
        xi_definite_integral = R(xi_integral(-1, eta) - xi_integral(-eta, eta))
        print(type(xi_definite_integral))
        eta_integral = R(xi_definite_integral.integral(eta))
        print(type(eta_integral))
        eta_definite_integral = eta_integral(xi, -1) - eta_integral(xi, 1)
        print(eta_definite_integral)
        return eta_definite_integral

    def inner_product(u, v):
        return integral_over_triangle(u * v) / 2

    def project(u, v):
        return inner_product(u, v) * v

    def gram_schmidt(v_list):
        e_list = []
        for v in v_list:
            u = v
            for e in e_list:
                u -= project(v, e)
            e = u / sa.sqrt(inner_product(u, u))
            e_list.append(e)
        return e_list

    max_degree = space_order - 1
    v_list = []
    for degree in range(max_degree + 1):
        for eta_degree in range(degree + 1):
            xi_degree = degree - eta_degree
            v_list.append(R((xi ** xi_degree) * (eta ** eta_degree)))
    print(v_list)
    e_list = gram_schmidt(v_list)

    coefficients = []
    for poly in e_list:
        poly_coefficients = []
        for degree in range(max_degree + 1):
            for eta_degree in range(degree + 1):
                xi_degree = degree - eta_degree
                poly_coefficients.append(poly.coefficient({xi: xi_degree, eta: eta_degree}))
        coefficients.append(poly_coefficients)

    return coefficients


def get_mass_matrix_2d_triangle(phi):
    return get_mass_matrix(
        phi, canonical_element.integrate_over_canonical_element_2d_triangle
    )

