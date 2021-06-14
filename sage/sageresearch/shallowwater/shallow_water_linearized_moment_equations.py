import sageresearch.shallowwater.shallow_water_moment_equations as swme

import sage.all as sa

#  Shallow Water Moment Equations
# See the paper
# Steady States and Well-balanced Schemes for Shallow Water Moment Equations with
# Topography by Koellermeier and Pimentel-Garcia


def get_swlme_equations_1d(num_moments):
    # Equations are of the form q_t + f(q)_x + G(q) q_x = 0
    # Quasilinear Form -> q_t + A q_x = 0
    num_eqns = num_moments + 2
    p = swme.get_primitive_variables_1d(num_moments)
    h = p[0]
    u = p[1]
    alpha = p[2:]

    list_ = swme.get_misc_variables()
    g = list_[0]
    e_z = list_[3]

    # Flux function
    f_list = [h * u]
    f_list.append(h * u * u + e_z * g * h * h / sa.Integer(2))
    if num_moments >= 1:
        f_list[1] += sum([sa.Integer(1) / sa.Integer(2 * i + 3) * h * alpha[i] * alpha[i] for i in range(num_moments)])
    for i in range(num_moments):
        f_list.append(sa.Integer(2) * h * u * alpha[i])

    f = sa.vector(f_list)

    # Nonconservative matrix
    G = sa.matrix(sa.SR, num_eqns)
    for i in range(num_moments):
        G[2 + i, 2 + i] = -u

    tuple_ = swme.get_substitution_dictionaries_1d(num_moments)
    p_to_q = tuple_[0]
    q_to_p = tuple_[1]
    q = swme.get_conserved_variables_1d(num_moments)

    # Quasilinear Matrix
    flux_jacobian = sa.jacobian(f.subs(p_to_q), q).subs(q_to_p)
    A = flux_jacobian + G

    return (f, G, A)


def get_swlme_eigenvalues_1d(num_moments):
    pass


def get_swlme_equations_2d(num_moments):
    num_eqns = 2 * num_moments + 3
    misc_var = swme.get_misc_variables()
    g = misc_var[0]
    e_z = misc_var[3]

    p = swme.get_primitive_variables_2d(num_moments)
    h = p[0]
    u = p[1]
    v = p[2]
    alpha = p[3::2]
    beta = p[4::2]

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

    for i in range(num_moments):
        f_x_list.append(
            sa.Integer(2) * h * u * alpha[i]
        )
        f_x_list.append(
            h * u * beta[i]
            + h * v * alpha[i]
        )

        f_y_list.append(
            h * u * beta[i]
            + h * v * alpha[i]
        )
        f_y_list.append(
            sa.Integer(2) * h * v * beta[i]
        )

    G_x_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]
    G_y_lists = [[sa.Integer(0) for i in range(num_eqns)] for j in range(num_eqns)]

    for i in range(num_moments):
        a_i_eqn = 3 + 2 * i
        b_i_eqn = 4 + 2 * i

        # (h alpha_i)_t + ... = u D_i + ...
        # u D_i = u (h alpha_i)_x + u (h beta_i)_y
        G_x_lists[a_i_eqn][a_i_eqn] = -u
        G_y_lists[a_i_eqn][b_i_eqn] = -u

        # (h beta_i)_t + ... = v D_i + ...
        # v D_i = v (h alpha_i)_x + v (h beta_i)_y
        G_x_lists[b_i_eqn][a_i_eqn] = -v
        G_y_lists[b_i_eqn][b_i_eqn] = -v

    f_x = sa.vector(f_x_list)
    f_y = sa.vector(f_y_list)
    G_x = sa.matrix(G_x_lists)
    G_y = sa.matrix(G_y_lists)

    tuple_ = swme.get_substitution_dictionaries_2d(num_moments)
    p_to_q = tuple_[0]
    q_to_p = tuple_[1]
    q = swme.get_conserved_variables_2d(num_moments)

    # flux_jacobians
    f_x_j = sa.jacobian(f_x.subs(p_to_q), q).subs(q_to_p)
    f_y_j = sa.jacobian(f_y.subs(p_to_q), q).subs(q_to_p)

    A_x = f_x_j + G_x
    A_y = f_y_j + G_y

    return (f_x, f_y, G_x, G_y, A_x, A_y)


def get_swlme_eigenvalues_2d(num_moments):
    pass
