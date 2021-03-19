# Hyperbolic Shallow Water Moment Equations
# See the paper
# Analysis and Numerical Simulation of Hyperbolic Shallow Water Moment Equations
# by Koellermeier and Rominger

import sageresearch.utils.symbolic_vector_matrix as symbolic_vector_matrix
import sageresearch.shallowwater.shallow_water_moment_equations as swme

import sage.all as sa


def get_hswme_equations(num_moments):
    num_eqns = num_moments + 2
    p = swme.get_primitive_variables_1d(num_moments)
    h = p[0]
    u = p[1]
    alpha = p[2:]

    g = sa.SR.symbol("g", domain="positive")
    e_z = sa.SR.symbol("e_z", domain="positive")

    # Quasilinear matrix
    A = sa.matrix(sa.SR, num_eqns)
    A[0, 1] = sa.Integer(1)
    A[1, 0] = e_z * g * h - u * u
    A[1, 1] = sa.Integer(2) * u
    if num_moments >= 1:
        A[1, 0] += -sa.Integer(1) / sa.Integer(3) * alpha[0] * alpha[0]
        A[1, 2] = sa.Integer(2) / sa.Integer(3) * alpha[0]
        A[2, 0] = -sa.Integer(2) * u * alpha[0]
        A[2, 1] = sa.Integer(2) * alpha[0]
        A[2, 2] = u
    if num_moments >= 2:
        A[3, 0] = -sa.Integer(2) / sa.Integer(3) * alpha[0] * alpha[0]
    for i in range(3, num_eqns):
        A[i - 1, i] = sa.Integer(i) / sa.Integer(2 * i - 1) * alpha[0]
        A[i, i] = u
        A[i, i - 1] = sa.Integer(i - 2) / sa.Integer(2 * i - 3) * alpha[0]

    return A


def get_hswme_eigenvalues(num_moments):
    # TODO: over six moments compute eigenvalues of A2 numerically
    p = swme.get_primitive_variables_1d(num_moments)
    h = p[0]
    u = p[1]
    misc_var = swme.get_misc_variables()
    g = misc_var[0]
    e_z = misc_var[3]
    if num_moments == 0:
        c = sa.sqrt(e_z * g * h)
        return [u + c, u - c]

    alpha_1 = p[2]
    c = sa.sqrt(e_z * g * h + alpha_1 * alpha_1)
    eigenvalues = [u + c, u - c]

    A_2 = sa.matrix(sa.SR, num_moments)
    for i in range(num_moments - 1):
        A_2[i, i + 1] = sa.Integer(i + 3) / sa.Integer(2 * i + 5)
        A_2[i + 1, i] = sa.Integer(i + 1) / sa.Integer(2 * i + 3)

    eig = A_2.eigenvalues()
    eigenvalues += [u + e.full_simplify() * alpha_1 for e in eig]

    return eigenvalues


def get_beta_hswme_equations(num_moments, beta=None):
    # num_moments is integer number of moments in vertical velocity profile
    # beta is list of symbolic expressions to add to final row
    # should be num_moments + 1 long
    if beta is None:
        # N - num_moments
        # default is beta_{N+1} = (N^2 - N)/(2N^2 + N - 1) alpha_1
        # beta_{0:N+1} = 0
        beta = [0 for i in range(num_moments + 1)]
        if num_moments >= 2:
            p = swme.get_primitive_variables_1d(num_moments)
            alpha_1 = p[2]
            beta[num_moments] = (
                sa.Integer(num_moments * num_moments - num_moments)
                / sa.Integer(2 * num_moments * num_moments + num_moments - 1)
                * alpha_1
            )

    A = get_hswme_equations(num_moments)
    for i in range(num_moments + 1):
        A[-1, i] += beta[i]

    return A


def get_beta_hswme_eigenvalues(num_moments, beta=None):
    pass
