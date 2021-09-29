from sageresearch.utils import legendre_polynomials

import sage.all as sa

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
    pass


def get_gauss_legendre_nodal_basis_1d(space_order):
    pass


def get_gauss_lobatto_nodal_basis_1d(space_order):
    pass


def get_legendre_basis_1d(space_order, inner_product_constant=one / two):
    return legendre_polynomials.get_legendre_polynomials(
        space_order, inner_product_constant
    )


def get_legendre_basis_2d_rectangular(space_order, inner_product_constant=one / four):
    phi = legendre_polynomials.get_legendre_polynomials(
        space_order, inner_product_constant
    )
    return phi


def get_mass_matrix_1d(basis):
    pass


def get_modal_basis_2d_triangle(space_order, inner_product_constant=one / two):
    pass


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

    coeffs = sa.matrix(15, 15)
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

    coeffs = coeffs / sa.sqrt(two * inner_product_constant)
    num_basis_cpts = int(space_order * (space_order + 1) / 2)
    return coeffs[:num_basis_cpts, :num_basis_cpts]
