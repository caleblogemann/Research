# Shallow Water Moment Equations in Spherical Coordinates
# See notes SphericalShallowWaterMomentsDerivation

import sageresearch.utils.symbolic_vector_matrix as svm
import sageresearch.utils.legendre_polynomials as legendre_polynomials
import sage.all as sa


def get_conserved_variables(num_moments):
    # get q
    num_eqns = 2 * num_moments + 3
    q = svm.get_vector_variable("q", num_eqns, "real")
    return q


def get_primitive_variables(num_moments):
    h = sa.SR.symbol("h", domain="positive")
    u_theta = sa.SR.symbol("u_theta", domain="real")
    u_phi = sa.SR.symbol("u_phi", domain="real")
    list_ = [h, u_theta, u_phi]
    if num_moments > 0:
        # add 1 to allow for 1 indexing
        alpha = svm.get_vector_variable("alpha", num_moments + 1, "real")
        beta = svm.get_vector_variable("beta", num_moments + 1, "real")

    for i in range(num_moments):
        list_.append(alpha[i + 1])
        list_.append(beta[i + 1])

    p = sa.vector(list_)
    return p


def get_substitution_dictionaries(num_moments):
    # get p_to_q, and q_to_p
    num_eqns = 2 * num_moments + 3
    q = get_conserved_variables(num_moments)
    p = get_primitive_variables(num_moments)
    h = p[0]
    u_theta = p[1]
    u_phi = p[2]
    p_to_q = {
        h: q[0],
        u_theta: q[1] / q[0],
        u_phi: q[2] / q[0],
    }
    q_to_p = {
        q[0]: h,
        q[1]: h * u_theta,
        q[2]: h * u_phi,
    }
    for i in range(num_moments):
        p_to_q[p[3 + 2 * i]] = q[3 + 2 * i] / q[0]
        p_to_q[p[4 + 2 * i]] = q[4 + 2 * i] / q[0]
        q_to_p[q[3 + 2 * i]] = h * p[3 + 2 * i]
        q_to_p[q[4 + 2 * i]] = h * p[4 + 2 * i]

    return (p_to_q, q_to_p)


def get_misc_variables():
    # return (g, e_x, e_y, e_z, nu, lambda_, h_b, t, x, y, z)
    t = sa.SR.symbol("t", domain="real")
    theta = sa.SR.symbol("theta", domain="real")
    phi = sa.SR.symbol("phi", domain="real")
    r = sa.SR.symbol("r", domain="positive")
    r_0 = sa.SR.symbol("r_0", domain="positive")
    g = sa.SR.symbol("g", domain="positive")
    e_theta = sa.SR.symbol("e_theta", domain="positive")
    e_phi = sa.SR.symbol("e_phi", domain="positive")
    e_r = sa.SR.symbol("e_r", domain="positive")
    nu = sa.SR.symbol("nu", domain="positive")
    lambda_ = sa.SR.symbol("lambda_", domain="positive")
    h_b = sa.function("h_b", nargs=2)(x, y)

    return sa.vector([g, e_theta, e_phi, e_r, nu, lambda_, h_b, t, theta, phi, r, r_0])


def get_shallow_water_moment_equations(num_moments):
    num_eqns = 3 + 2 * num_moments
    misc_var = get_misc_variables()
    g = misc_var[0]
    e_theta = misc_var[1]
    e_phi = misc_var[2]
    e_r = misc_var[3]
    t = misc_var[7]
    theta = misc_var[8]
    phi = misc_var[9]
    r_0 = misc_var[11]

    p = get_primitive_variables_2d(num_moments)
    h = p[0]
    u_theta = p[1]
    u_phi = p[2]
    alpha = p[3::2]
    beta = p[4::2]

    psi = get_velocity_basis_functions(num_moments)

    A = [
        [
            [
                sa.Integer(2 * i + 3)
                * (phi[i + 1] * phi[j + 1] * phi[k + 1]).integrate(
                    x, sa.Integer(0), sa.Integer(1)
                )
                for k in range(num_moments)
            ]
            for j in range(num_moments)
        ]
        for i in range(num_moments)
    ]
    B = [
        [
            [
                sa.Integer(2 * i + 3)
                * (
                    phi[i + 1].derivative(x)
                    * (phi[j + 1](x=z)).integrate(z, sa.Integer(0), x)
                    * phi[k + 1]
                ).integrate(x, 0, 1)
                for k in range(num_moments)
            ]
            for j in range(num_moments)
        ]
        for i in range(num_moments)
    ]

    sum_a2 = sum(
        [alpha[j] * alpha[j] / sa.Integer(2 * j + 3) for j in range(num_moments)]
    )
    sum_b2 = sum(
        [beta[j] * beta[j] / sa.Integer(2 * j + 3) for j in range(num_moments)]
    )
    sum_ab = sum(
        [alpha[j] * beta[j] / sa.Integer(2 * j + 3) for j in range(num_moments)]
    )

    r_sin_phi = r_0 * sa.sin(phi)

    f_theta_list = [
        h * u_theta / r_sin_phi,
        (
            h * u_theta * u_theta
            + h * sum_a2
            + sa.Integer(1) / sa.Integer(2) * g * e_r * h * h
        )
        / r_sin_phi,
        (h * u_theta * u_phi + h * sum_ab) / r_sin_phi,
    ]
    f_phi_list = [
        h * u_phi / r_0,
        (h * u_theta * u_phi + h * sum_ab) / r_0,
        (
            h * u_phi * u_phi
            + h * sum_b2
            + sa.Integer(1) / sa.Integer(2) * g * e_r * h * h
        )
        / r_0,
    ]
    s_list = [
        h * u_phi * sa.cot(phi) / r_0,
        -sa.Integer(2) * h * sa.cot(phi) / r * (u_theta * u_phi + sum_ab)
        - g * h * e_r * h_b.diff(theta) / r_sin_phi
        + g * h * e_theta,
        -h * sa.cot(phi) / r_0 * (u_phi * u_phi - u_theta * u_theta + sum_b2 - sum_a2)
        - g * h * e_r / r_0 * h_b.diff(phi)
        + g * h * e_phi,
    ]

    G_phi

    for i in range(num_moments):
        # h alpha_i equation
        f_theta_list.append()
        f_phi_list.append()
        s_list.append()

        # h beta_i equation
        f_theta_list.append()
        f_phi_list.append()
        s_list.append()


def get_velocity_basis_functions(num_moments):
    # velocity is u \phi_0 + \sum{j=1}{num_moments}{alpha_j \phi_j(\zeta)}
    # return functions phi
    phi = legendre_polynomials.get_legendre_polynomials_fixed_lower_endpoint(
        num_moments + 1, sa.Integer(1), sa.Integer(0), sa.Integer(1)
    )
    return phi
