# Shallow Water Moment Equations
# See the paper, Shallow Water Moment Cascades
# by Kowalski and Torrilhon

import sageresearch.utils.symbolic_vector_matrix as svm
import sageresearch.utils.legendre_polynomials as legendre_polynomials

import sage.all as sa


def get_conserved_variables_1d(num_moments):
    # get q
    num_eqns = num_moments + 2
    q = svm.get_vector_variable("q", num_eqns, "real")
    return q


def get_primitive_variables_1d(num_moments):
    # get p
    h = sa.SR.symbol("h", domain="positive")
    u = sa.SR.symbol("u", domain="real")
    list_ = [h, u]
    if num_moments > 0:
        # add 1 to allow for 1 indexing
        alpha = svm.get_vector_variable("alpha", num_moments + 1, "real")
        list_ += alpha[1:].list()

    p = sa.vector(list_)
    return p


def get_substitution_dictionaries_1d(num_moments):
    # get p_to_q, and q_to_p
    num_eqns = num_moments + 2
    q = get_conserved_variables_1d(num_moments)
    p = get_primitive_variables_1d(num_moments)
    h = p[0]
    u = p[1]
    p_to_q = {
        h: q[0],
        u: q[1] / q[0],
    }
    q_to_p = {
        q[0]: h,
        q[1]: h * u,
    }
    for i in range(num_moments):
        p_to_q[p[2 + i]] = q[2 + i] / q[0]
        q_to_p[q[2 + i]] = h * p[2 + i]

    return (p_to_q, q_to_p)
    pass


def get_conserved_variables_2d(num_moments):
    num_eqns = 2 * num_moments + 3
    q = svm.get_vector_variable("q", num_eqns, "real")
    return q


def get_primitive_variables_2d(num_moments):
    h = sa.SR.symbol("h", domain="positive")
    u = sa.SR.symbol("u", domain="real")
    v = sa.SR.symbol("v", domain="real")
    if num_moments > 0:
        # add 1 to allow for 1 indexing
        alpha = svm.get_vector_variable("alpha", num_moments + 1, "real")
        beta = svm.get_vector_variable("beta", num_moments + 1, "real")

    list_ = [h, u, v]
    for i in range(num_moments):
        list_.append(alpha[i + 1])
        list_.append(beta[i + 1])
    return sa.vector(list_)


def get_substitution_dictionaries_2d(num_moments):
    # get p_to_q, and q_to_p
    num_eqns = 2 * num_moments + 3
    q = get_conserved_variables_2d(num_moments)
    p = get_primitive_variables_2d(num_moments)
    h = p[0]
    u = p[1]
    v = p[2]
    p_to_q = {
        h: q[0],
        u: q[1] / q[0],
        v: q[2] / q[0],
    }
    q_to_p = {
        q[0]: h,
        q[1]: h * u,
        q[2]: h * v,
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
    x = sa.SR.symbol("x", domain="real")
    y = sa.SR.symbol("y", domain="real")
    z = sa.SR.symbol("z", domain="real")
    g = sa.SR.symbol("g", domain="positive")
    e_x = sa.SR.symbol("e_x", domain="positive")
    e_y = sa.SR.symbol("e_y", domain="positive")
    e_z = sa.SR.symbol("e_z", domain="positive")
    nu = sa.SR.symbol("nu", domain="positive")
    lambda_ = sa.SR.symbol("lambda_", domain="positive")
    h_b = sa.function("h_b", nargs=2)(x, y)

    return sa.vector([g, e_x, e_y, e_z, nu, lambda_, h_b, t, x, y, z])

# TODO: potentially add is_functions
def get_generalized_shallow_water_equations_1d(num_moments):
    # get f_p, G_p, s_p, A_p

    p = get_primitive_variables_2d(num_moments)
    # substitution dict to set v, and all betas to zero
    dict_ = {p[2]: 0}
    for i in range(num_moments):
        dict_[p[4 + 2 * i]] = 0

    # get 2d equations (f_x_p, f_y_p, G_x_p, G_y_p, s_p)
    tuple_ = get_generalized_shallow_water_equations_2d(num_moments)

    # set to zero in x variables
    f_p = tuple_[0].subs(dict_)
    G_p = tuple_[2].subs(dict_)
    s_p = tuple_[4].subs(dict_)
    A_p = tuple_[5].subs(dict_)

    # delete rows and and colums related to v and betas
    rows_delete = [2 * i + 2 for i in range(num_moments + 1)]
    G_p = G_p.delete_rows(rows_delete)
    G_p = G_p.delete_columns(rows_delete)

    A_p = A_p.delete_rows(rows_delete)
    A_p = A_p.delete_columns(rows_delete)

    rows_keep = [2 * i + 1 for i in range(num_moments + 1)]
    rows_keep.insert(0, 0)
    f_p = sa.vector([f_p[i] for i in rows_keep])
    s_p = sa.vector([s_p[i] for i in rows_keep])

    return (f_p, G_p, s_p, A_p)


def get_generalized_shallow_water_equations_2d(num_moments):
    # get f_x_p, f_y_p, G_x_p, G_y_p, s_p, A_x_p, A_y_p
    num_eqns = 3 + 2 * num_moments
    misc_var = get_misc_variables()
    g = misc_var[0]
    e_x = misc_var[1]
    e_y = misc_var[2]
    e_z = misc_var[3]
    nu = misc_var[4]
    lambda_ = misc_var[5]
    h_b = misc_var[6]
    t = misc_var[7]
    x = misc_var[8]
    y = misc_var[9]
    z = misc_var[10]

    p = get_primitive_variables_2d(num_moments)
    h = p[0]
    u = p[1]
    v = p[2]
    alpha = p[3::2]
    beta = p[4::2]

    phi = legendre_polynomials.get_legendre_polynomials_fixed_lower_endpoint(
        num_moments + 1, sa.Integer(1), sa.Integer(0), sa.Integer(1)
    )
    A = [
        [
            [
                sa.Integer(2 * i + 3) * (phi[i + 1] * phi[j + 1] * phi[k + 1]).integrate(x, sa.Integer(0), sa.Integer(1))
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
    C = [
        [
            (phi[i + 1].derivative(x) * phi[j + 1].derivative(x)).integrate(x, sa.Integer(0), sa.Integer(1))
            for j in range(num_moments)
        ]
        for i in range(num_moments)
    ]
    D = [
        (h * alpha[i]).derivative(x) + (h * beta[i]).derivative(y)
        for i in range(num_moments - 1)
    ]

    # fluxes
    f_x_list = [
        h * u,
        h * u * u
        + h * sum([alpha[j] * alpha[j] / sa.Integer(2 * j + 3) for j in range(num_moments)])
        + sa.Integer(1) / sa.Integer(2) * g * e_z * h * h,
        h * u * v
        + h * sum([alpha[j] * beta[j] / sa.Integer(2 * j + 3) for j in range(num_moments)]),
    ]
    f_y_list = [
        h * v,
        h * u * v
        + h * sum([alpha[j] * beta[j] / sa.Integer(2 * j + 3) for j in range(num_moments)]),
        h * v * v
        + h * sum([beta[j] * beta[j] / sa.Integer(2 * j + 3) for j in range(num_moments)])
        + sa.Integer(1) / sa.Integer(2) * g * e_z * h * h,
    ]

    # source term
    s_list = [
        sa.Integer(0),
        -nu / lambda_ * (u + sum([alpha[j] for j in range(num_moments)]))
        + h * g * (e_x - e_z * h_b.derivative(x)),
        -nu / lambda_ * (v + sum([beta[j] for j in range(num_moments)]))
        + h * g * (e_y - e_z * h_b.derivative(y)),
    ]

    for i in range(num_moments):
        f_x_list.append(
            sa.Integer(2) * h * u * alpha[i]
            + h
            * sum(
                [
                    sum([A[i][j][k] * alpha[j] * alpha[k] for k in range(num_moments)])
                    for j in range(num_moments)
                ]
            )
        )
        f_x_list.append(
            h * u * beta[i]
            + h * v * alpha[i]
            + h
            * sum(
                [
                    sum([A[i][j][k] * alpha[j] * beta[k] for k in range(num_moments)])
                    for j in range(num_moments)
                ]
            )
        )

        f_y_list.append(
            h * u * beta[i]
            + h * v * alpha[i]
            + h
            * sum(
                [
                    sum([A[i][j][k] * alpha[j] * beta[j] for k in range(num_moments)])
                    for j in range(num_moments)
                ]
            )
        )
        f_y_list.append(
            sa.Integer(2) * h * v * beta[i]
            + h
            * sum(
                [
                    sum([A[i][j][k] * beta[j] * beta[k] for k in range(num_moments)])
                    for j in range(num_moments)
                ]
            )
        )

        s_list.append(
            -sa.Integer(2 * i + 3)
            * nu
            / lambda_
            * (
                u
                + sum(
                    [(sa.Integer(1) + lambda_ / h * C[i][j]) * alpha[j] for j in range(num_moments)]
                )
            )
        )
        s_list.append(
            -sa.Integer(2 * i + 3)
            * nu
            / lambda_
            * (
                v
                + sum(
                    [(sa.Integer(1) + lambda_ / h * C[i][j]) * beta[j] for j in range(num_moments)]
                )
            )
        )

    # nonconservative product matrices
    G_x_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]
    G_y_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]

    for i in range(num_moments):
        a_i_eqn = 3 + 2 * i
        b_i_eqn = 4 + 2 * i
        for j in range(num_moments):
            # (h alpha_j)_x term, column index
            a_j_eqn = 3 + 2 * j
            # (h beta_j)_y term, column index
            b_j_eqn = 4 + 2 * j

            # h alpha[i] equation
            # (h alpha_i)_t + ... = ... - sum{j=1}{N}{D_j \sum{k=1}{N}{B_ijk alpha_k}} + ...
            # (h alpha_j)_x sum{k=1}{N}{-B_ijk alpha_k}
            G_x_lists[a_i_eqn][a_j_eqn] = sum(
                [-B[i][j][k] * alpha[k] for k in range(num_moments)]
            )
            # (h beta_j)_y sum{k=1}{N}{-B_ijk alpha_k}
            G_y_lists[a_i_eqn][b_j_eqn] = sum(
                [-B[i][j][k] * alpha[k] for k in range(num_moments)]
            )

            # (h beta_i)_t + ... = ... - sum{j=1}{N}{D_j \sum{k=1}{N}{B_ijk beta_k}} + ...
            # (h alpha_j)_x sum{k=1}{N}{-B_ijk beta_k}
            G_x_lists[b_i_eqn][a_j_eqn] = sum(
                [-B[i][j][k] * beta[k] for k in range(num_moments)]
            )
            # (h beta_j)_y sum{k=1}{N}{-B_ijk beta_k}
            G_y_lists[b_i_eqn][b_j_eqn] = sum(
                [-B[i][j][k] * beta[k] for k in range(num_moments)]
            )

        # (h alpha_i)_t + ... = u D_i + ...
        # u D_i = u (h alpha_i)_x + u (h beta_i)_y
        G_x_lists[a_i_eqn][a_i_eqn] += u
        G_y_lists[a_i_eqn][b_i_eqn] += u

        # (h beta_i)_t + ... = v D_i + ...
        # v D_i = v (h alpha_i)_x + v (h beta_i)_y
        G_x_lists[b_i_eqn][a_i_eqn] += v
        G_y_lists[b_i_eqn][b_i_eqn] += v

    f_x_p = sa.vector(f_x_list)
    f_y_p = sa.vector(f_y_list)
    G_x_p = sa.matrix(G_x_lists)
    G_y_p = sa.matrix(G_y_lists)
    s_p = sa.vector(s_list)

    tuple_ = get_substitution_dictionaries_2d(num_moments)
    p_to_q = tuple_[0]
    q_to_p = tuple_[1]
    q = get_conserved_variables_2d(num_moments)

    f_x_q = f_x_p.subs(p_to_q)
    f_y_q = f_y_p.subs(p_to_q)

    # flux_jacobians
    f_x_q_j = sa.jacobian(f_x_q, q)
    f_y_q_j = sa.jacobian(f_y_q, q)

    # Nonconserved matrix
    G_x_q = G_x_p.subs(p_to_q)
    G_y_q = G_y_p.subs(p_to_q)

    # Quasilinear Matrices
    A_x_q = f_x_q_j - G_x_q
    A_y_q = f_y_q_j - G_y_q
    # should also be equivalent to f_x_p_j @ p_q_j - G_x_p
    A_x_p = A_x_q.subs(q_to_p)
    # should also be equivalent to f_y_p_j @ p_q_j - G_y_p
    A_y_p = A_y_q.subs(q_to_p)

    return (f_x_p, f_y_p, G_x_p, G_y_p, s_p, A_x_p, A_y_p)
