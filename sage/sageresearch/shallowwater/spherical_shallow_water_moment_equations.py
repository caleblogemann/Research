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
    h_b = sa.function("h_b", nargs=2)(theta, phi)

    return sa.vector([g, e_theta, e_phi, e_r, nu, lambda_, h_b, t, theta, phi, r, r_0])


def get_shallow_water_moment_equations(num_moments):
    num_eqns = 3 + 2 * num_moments
    misc_var = get_misc_variables()
    g = misc_var[0]
    e_theta = misc_var[1]
    e_phi = misc_var[2]
    e_r = misc_var[3]
    h_b = misc_var[6]
    t = misc_var[7]
    theta = misc_var[8]
    phi = misc_var[9]
    r = misc_var[10]
    r_0 = misc_var[11]

    p = get_primitive_variables(num_moments)
    h = p[0]
    u_theta = p[1]
    u_phi = p[2]
    alpha = p[3::2]
    beta = p[4::2]

    tuple_ = get_velocity_basis_functions(num_moments)
    psi = tuple_[0]
    x = tuple_[1]

    A = [
        [
            [
                sa.Integer(2 * i + 3)
                * (psi[i + 1] * psi[j + 1] * psi[k + 1]).integrate(
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
                    psi[i + 1].derivative(x)
                    * (psi[k + 1](x=r)).integrate(r, sa.Integer(0), x)
                    * psi[j + 1]
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
    h_cot_phi_r = h * sa.cot(phi) / r_0

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
        -h_cot_phi_r * u_phi,
        -sa.Integer(2) * h_cot_phi_r * (u_theta * u_phi + sum_ab)
        - g * h * e_r * h_b.derivative(theta) / r_sin_phi
        + g * h * e_theta,
        -h_cot_phi_r * (u_phi * u_phi - u_theta * u_theta + sum_b2 - sum_a2)
        - g * h * e_r / r_0 * h_b.derivative(phi)
        + g * h * e_phi,
    ]

    G_theta_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]
    G_phi_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]

    for i in range(num_moments):
        sum_Aaa = sum(
            [
                sum([alpha[j] * alpha[k] * A[i][j][k] for k in range(num_moments)])
                for j in range(num_moments)
            ]
        )
        sum_Aab = sum(
            [
                sum([alpha[j] * beta[k] * A[i][j][k] for k in range(num_moments)])
                for j in range(num_moments)
            ]
        )
        sum_Abb = sum(
            [
                sum([A[i][j][k] * beta[j] * beta[k] for k in range(num_moments)])
                for j in range(num_moments)
            ]
        )
        sum_Bab = sum(
            [
                sum([B[i][j][k] * alpha[j] * beta[k] for k in range(num_moments)])
                for j in range(num_moments)
            ]
        )
        sum_Bbb = sum(
            [
                sum([B[i][j][k] * beta[j] * beta[k] for k in range(num_moments)])
                for j in range(num_moments)
            ]
        )
        # index for alpha_i equation and beta_i equation
        a_i_eqn = 3 + 2 * i
        b_i_eqn = 4 + 2 * i

        # h alpha_i equation
        f_theta_list.append(
            (sa.Integer(2) * h * u_theta * alpha[i] + h * sum_Aaa) / r_sin_phi
        )
        f_phi_list.append(
            (h * u_theta * beta[i] + h * u_phi * alpha[i] + h * sum_Aab) / r_0
        )
        s_list.append(
            -h_cot_phi_r
            * (
                u_theta * beta[i]
                + sa.Integer(2) * u_phi * alpha[i]
                + sa.Integer(2) * sum_Aab
                + sum_Bab
            )
        )
        G_theta_lists[a_i_eqn][a_i_eqn] += -u_theta / r_sin_phi
        G_phi_lists[a_i_eqn][b_i_eqn] += -u_theta / r_0

        # h beta_i equation
        f_theta_list.append(
            (h * u_theta * beta[i] + h * u_phi * alpha[i] + h * sum_Aab) / r_sin_phi
        )
        f_phi_list.append((sa.Integer(2) * h * u_phi * beta[i] + h * sum_Abb) / r_0)
        s_list.append(
            h_cot_phi_r
            * (
                sa.Integer(2) * u_theta * alpha[i]
                - u_phi * beta[i]
                + sum_Aaa
                - sum_Abb
                - sum_Bbb
            )
        )
        G_theta_lists[b_i_eqn][a_i_eqn] += -u_phi / r_sin_phi
        G_phi_lists[b_i_eqn][b_i_eqn] += -u_phi / r_0

        for k in range(num_moments):
            a_k_eqn = 3 + 2 * k
            b_k_eqn = 4 + 2 * k

            sum_Ba = sum([alpha[j] * B[i][j][k] for j in range(num_moments)])
            sum_Bb = sum([beta[j] * B[i][j][k] for j in range(num_moments)])

            # h alpha_i equation
            G_theta_lists[a_i_eqn][a_k_eqn] += sum_Ba / r_sin_phi
            G_phi_lists[a_i_eqn][b_k_eqn] += sum_Ba / r_0

            # h beta_i equation
            G_theta_lists[b_i_eqn][a_k_eqn] += sum_Bb / r_sin_phi
            G_phi_lists[b_i_eqn][b_k_eqn] += sum_Bb / r_0

    f_theta_p = sa.vector(f_theta_list)
    f_phi_p = sa.vector(f_phi_list)
    G_theta_p = sa.matrix(G_theta_lists)
    G_phi_p = sa.matrix(G_phi_lists)
    s_p = sa.vector(s_list)

    tuple_ = get_substitution_dictionaries(num_moments)
    p_to_q = tuple_[0]
    q_to_p = tuple_[1]
    q = get_conserved_variables(num_moments)

    f_theta_q = f_theta_p.subs(p_to_q)
    f_phi_q = f_phi_p.subs(p_to_q)

    #flux jacobians
    f_theta_q_j = sa.jacobian(f_theta_q, q)
    f_phi_q_j = sa.jacobian(f_phi_q, q)

    G_theta_q = G_theta_p.subs(p_to_q)
    G_phi_q = G_phi_p.subs(p_to_q)

    # Quasilinear Matrix
    A_theta_q = f_theta_q_j + G_theta_q
    A_phi_q = f_phi_q_j + G_phi_q
    A_theta_p = A_theta_q.subs(q_to_p)
    A_phi_p = A_phi_q.subs(q_to_p)

    return (f_theta_p, f_phi_p, G_theta_p, G_phi_p, s_p, A_theta_p, A_phi_p)


def get_velocity_basis_functions(num_moments):
    # velocity is u \phi_0 + \sum{j=1}{num_moments}{alpha_j \phi_j(\zeta)}
    # return functions phi
    phi = legendre_polynomials.get_legendre_polynomials_fixed_lower_endpoint(
        num_moments + 1, sa.Integer(1), sa.Integer(0), sa.Integer(1)
    )
    x = sa.var("x")
    return (phi, x)
