from sageresearch.utils import legendre_polynomials

import sage.all as sa


def get_nodal_basis_1d(nodes):
    pass


def get_gauss_legendre_nodal_basis_1d(space_order):
    pass


def get_gauss_lobatto_nodal_basis_1d(space_order):
    pass


def get_legendre_basis_1d(space_order, inner_product_constant=0.5):
    return legendre_polynomials.get_legendre_polynomials(space_order, inner_product_constant)


def get_legendre_basis_2d_rectangular(space_order, inner_product_constant=0.25):
    phi = legendre_polynomials.get_legendre_polynomials(space_order, inner_product_constant)
    return phi


def get_mass_matrix_1d(basis):
    pass
