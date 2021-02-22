load("../Sage/SymbolicVectorMatrix.sage")
load("../Sage/LegendrePolynomials.sage")


def set_generalized_shallow_water_variables_1d(num_moments):
    g = var("g", domain="positive")
    # primitive variables
    h = var("h", domain="positive")
    u = var("u")
    if num_moments > 0:
        alpha = get_vector_variable("alpha", num_moments)

    list_ = [h, u]
    for i in range(num_moments):
        list_.append(alpha[i])

    global primitive, conserved
    primitive = vector(list_)
    conserved = get_vector_variable("q", 2 + num_moments)

    global p, q
    p = primitive
    q = conserved

    # substitution dictionaries
    global p_to_q, q_to_p
    p_to_q = {
        h: q[0],
        u: q[1] / q[0],
    }
    q_to_p = {
        q[0]: h,
        q[1]: h * u,
    }
    for i in range(num_moments):
        p_to_q[alpha[i]] = q[2 + i] / q[0]
        q_to_p[q[2 + i]] = h * alpha[i]

    # transformation symbolic expressions
    global p_q, q_p
    p_q = vector(p_to_q.values())
    q_p = vector(q_to_p.values())

    # transformation jacobians
    global p_q_j, q_p_j
    p_q_j = jacobian(p_q, q)
    q_p_j = jacobian(q_p, p)

    # get equations, get f_p, G_p, s_p
    get_generalized_shallow_water_equations_1d(num_moments)

    # conserved and transform fluxes
    global f_q
    f_q = f_p.subs(p_to_q)

    # flux_jacobians
    global f_p_j, f_q_j
    f_p_j = jacobian(f_p, p)
    f_q_j = jacobian(f_q, q)

    # Nonconserved matrix
    global G_q
    G_q = G_p.subs(p_to_q)

    # Quasilinear Matrix
    global A_p, A_q
    A_q = f_q_j - G_q
    # should also be equivalent to f_p_j @ p_q_j - G_p
    A_p = A_q.subs(q_to_p)

    # TODO: needs to be rewritten/rechecked
    # transform variables
    # global transform
    # transform = get_vector_variable('z', num_moments + 2)

    # global z
    # z = transform

    # substitution dictionaries
    # global p_to_z, q_to_z, z_to_p, z_to_q
    # primitive_to_transform = {
        # h: z_1 ^ 2,
        # u: z_2 / z_1,
        # s: z_3 / z_1,
        # k: z_4 / z_1,
        # m: z_5 / z_1,
    # }
    # p_to_z = primitive_to_transform
    # conserved_to_transform = {
        # q_1: z_1 ^ 2,
        # q_2: z_1 * z_2,
        # q_3: z_1 * z_3,
        # q_4: z_1 * z_4,
        # q_5: z_1 * z_5,
    # }
    # q_to_z = conserved_to_transform
    # transform_to_primitive = {
        # z_1: sqrt(h),
        # z_2: sqrt(h) * u,
        # z_3: sqrt(h) * s,
        # z_4: sqrt(h) * k,
        # z_5: sqrt(h) * m,
    # }
    # z_to_p = transform_to_primitive
    # transform_to_conserved = {
        # z_1: sqrt(q_1),
        # z_2: q_2 / sqrt(q_1),
        # z_3: q_3 / sqrt(q_1),
        # z_4: q_4 / sqrt(q_1),
        # z_5: q_5 / sqrt(q_1),
    # }
    # z_to_q = transform_to_conserved

    # transformation symbolic expressions
    # global p_z, q_z, z_p, z_q
    # p(z)
    # temp = vector([z_1 ^ 2, z_2 / z_1, z_3 / z_1, z_4 / z_1, z_5 / z_1])
    # primitive_from_transform = temp[0 : (num_moments + 2)]
    # p_z = primitive_from_transform

    # q(z)
    # temp = vector([z_1 ^ 2, z_1 * z_2, z_1 * z_3, z_1 * z_4, z_1 * z_4])
    # conserved_from_transform = temp[0 : (num_moments + 2)]
    # q_z = conserved_from_transform

    # z(p)
    # temp = vector([sqrt(h), sqrt(h) * u, sqrt(h) * s, sqrt(h) * k, sqrt(h) * m])
    # transform_from_primitive = temp[0 : (num_moments + 2)]
    # z_p = transform_from_primitive
    # z(q)
    # temp = vector(
        # [sqrt(q_1), q_2 / sqrt(q_1), q_3 / sqrt(q_1), q_4 / sqrt(q_1), q_5 / sqrt(q_1)]
    # )
    # transform_from_conserved = temp[0 : (num_moments + 2)]
    # z_q = transform_from_conserved

    # transformation jacobians
    # global  p_z_j, q_z_j, z_p_j, z_q_j
    # p_z_j = jacobian(p_z, z)
    # q_z_j = jacobian(q_z, z)
    # z_p_j = jacobian(z_p, p)
    # z_q_j = jacobian(z_q, q)

    # conserved and transform fluxes
    # global flux_transform, f_z
    # flux_transform = flux_primitive.subs(p_to_z)
    # f_z = flux_transform

    # flux_jacobians
    # global f_z_j
    # global flux_transform_jacobian
    # flux_transform_jacobian = jacobian(f_z, z)
    # f_z_j = flux_transform_jacobian

    # Nonconserved matrix
    # global Q_q, Q_z
    # Q_q = Q_p.subs(p_to_q)
    # Q_z = Q_p.subs(p_to_z)

    # Quasilinear Matrices
    # global A_z
    # A_z = A_q.subs(q_to_z)

    # Riemann Problem
    # global p_to_P_i, P_i, p_to_P_im1, P_im1
    # P_i = [H_i, U_i]
    # H_i = var("H_i", domain="positive")
    # U_i = var("U_i")
    # S_i = var("S_i")
    # K_i = var("K_i")
    # M_i = var("M_i")
    # p_to_P_i = {h: H_i, u: U_i, s: S_i, k: K_i, m: M_i}
    # P_i = p.subs(p_to_P_i)
    # P_im1 = [H_im1, U_im1]
    # H_im1 = var("H_im1", domain="positive")
    # U_im1 = var("U_im1")
    # S_im1 = var("S_im1")
    # K_im1 = var("K_im1")
    # M_im1 = var("M_im1")
    # p_to_P_im1 = {h: H_im1, u: U_im1, s: S_im1}
    # P_im1 = p.subs(p_to_P_im1)
    # global q_to_Q_i, Q_i, q_to_Q_im1, Q_im1
    # Q_i = [Q_1_i, Q_2_i]
    # Q_1_i = var("Q_1_i", domain="positive")
    # Q_2_i = var("Q_2_i")
    # Q_3_i = var("Q_3_i")
    # Q_4_i = var("Q_4_i")
    # Q_5_i = var("Q_5_i")
    # q_to_Q_i = {q_1: Q_1_i, q_2: Q_2_i, q_3: Q_3_i, q_4: Q_4_i, q_5: Q_5_i}
    # Q_i = q.subs(q_to_Q_i)
    # Q_im1 = [Q_1_im1, Q_2_im1]
    # Q_1_im1 = var("Q_1_im1", domain="positive")
    # Q_2_im1 = var("Q_2_im1")
    # Q_3_im1 = var("Q_3_im1")
    # Q_4_im1 = var("Q_4_im1")
    # Q_5_im1 = var("Q_5_im1")
    # q_to_Q_im1 = {q_1: Q_1_im1, q_2: Q_2_im1, q_3: Q_3_im1, q_4: Q_4_im1, q_5: Q_5_im1}
    # Q_im1 = q.subs(q_to_Q_im1)
    # global z_to_Z_i, Z_i, z_to_Z_im1, Z_im1
    # Z_i = [Z_1_i, Z_2_i]
    # Z_1_i = var("Z_1_i", domain="positive")
    # Z_2_i = var("Z_2_i")
    # Z_3_i = var("Z_3_i")
    # Z_4_i = var("Z_4_i")
    # Z_5_i = var("Z_5_i")
    # z_to_Z_i = {z_1: Z_1_i, z_2: Z_2_i, z_3: Z_3_i, z_4: Z_4_i, z_5: Z_5_i}
    # Z_i = z.subs(z_to_Z_i)
    # Z_im1 = [Z_1_im1, Z_2_im1]
    # Z_1_im1 = var("Z_1_im1", domain="positive")
    # Z_2_im1 = var("Z_2_im1")
    # Z_3_im1 = var("Z_3_im1")
    # Z_4_im1 = var("Z_4_im1")
    # Z_5_im1 = var("Z_5_im1")
    # z_to_Z_im1 = {z_1: Z_1_im1, z_2: Z_2_im1, z_3: Z_3_im1, z_4: Z_4_im1, z_5: Z_5_im1}
    # Z_im1 = z.subs(z_to_Z_im1)

    # global P_i_to_Q_i, P_i_to_Z_i
    # global Q_i_to_P_i, Q_i_to_Z_i
    # global Z_i_to_P_i, Z_i_to_Q_i
    # P_i_to_Q_i = {
        # H_i: Q_1_i,
        # U_i: Q_2_i / Q_1_i,
        # S_i: Q_3_i / Q_1_i,
        # K_i: Q_4_i / Q_1_i,
        # M_i: Q_5_i / Q_1_i,
    # }
    # P_i_to_Z_i = {
        # H_i: Z_1_i ^ 2,
        # U_i: Z_2_i / Z_1_i,
        # S_i: Z_3_i / Z_1_i,
        # K_i: Z_4_i / Z_1_i,
        # M_i: Z_5_i / Z_1_i,
    # }
    # Q_i_to_P_i = {
        # Q_1_i: H_i,
        # Q_2_i: H_i * U_i,
        # Q_3_i: H_i * S_i,
        # Q_4_i: H_i * K_i,
        # Q_5_i: H_i * M_i,
    # }
    # Q_i_to_Z_i = {
        # Q_1_i: Z_1_i ^ 2,
        # Q_2_i: Z_1_i * Z_2_i,
        # Q_3_i: Z_1_i * Z_3_i,
        # Q_4_i: Z_1_i * Z_4_i,
        # Q_5_i: Z_1_i * Z_5_i,
    # }
    # Z_i_to_Q_i = {
        # Z_1_i: sqrt(Q_1_i),
        # Z_2_i: Q_2_i / sqrt(Q_1_i),
        # Z_3_i: Q_3_i / sqrt(Q_1_i),
        # Z_4_i: Q_4_i / sqrt(Q_1_i),
        # Z_5_i: Q_5_i / sqrt(Q_1_i),
    # }
    # Z_i_to_P_i = {
        # Z_1_i: sqrt(H_i),
        # Z_2_i: sqrt(H_i) * U_i,
        # Z_3_i: sqrt(H_i) * S_i,
        # Z_4_i: sqrt(H_i) * K_i,
        # Z_5_i: sqrt(H_i) * M_i,
    # }

    # global P_im1_to_Q_im1, P_im1_to_Z_im1
    # global Q_im1_to_P_im1, Q_im1_to_Z_im1
    # global Z_im1_to_P_im1, Z_im1_to_Q_im1
    # P_im1_to_Q_im1 = {
        # H_im1: Q_1_im1,
        # U_im1: Q_2_im1 / Q_1_im1,
        # S_im1: Q_3_im1 / Q_1_im1,
        # K_im1: Q_4_im1 / Q_1_im1,
        # M_im1: Q_5_im1 / Q_1_im1,
    # }
    # P_im1_to_Z_im1 = {
        # H_im1: Z_1_im1 ^ 2,
        # U_im1: Z_2_im1 / Z_1_im1,
        # S_im1: Z_3_im1 / Z_1_im1,
        # K_im1: Z_4_im1 / Z_1_im1,
        # M_im1: Z_5_im1 / Z_1_im1,
    # }
    # Q_im1_to_P_im1 = {
        # Q_1_im1: H_im1,
        # Q_2_im1: H_im1 * U_im1,
        # Q_3_im1: H_im1 * S_im1,
        # Q_4_im1: H_im1 * K_im1,
        # Q_5_im1: H_im1 * M_im1,
    # }
    # Q_im1_to_Z_im1 = {
        # Q_1_im1: Z_1_im1 ^ 2,
        # Q_2_im1: Z_1_im1 * Z_2_im1,
        # Q_3_im1: Z_1_im1 * Z_3_im1,
        # Q_4_im1: Z_1_im1 * Z_4_im1,
        # Q_5_im1: Z_1_im1 * Z_5_im1,
    # }
    # Z_im1_to_Q_im1 = {
        # Z_1_im1: sqrt(Q_1_im1),
        # Z_2_im1: Q_2_im1 / sqrt(Q_1_im1),
        # Z_3_im1: Q_3_im1 / sqrt(Q_1_im1),
        # Z_4_im1: Q_4_im1 / sqrt(Q_1_im1),
        # Z_5_im1: Q_5_im1 / sqrt(Q_1_im1),
    # }
    # Z_im1_to_P_im1 = {
        # Z_1_im1: sqrt(H_im1),
        # Z_2_im1: sqrt(H_im1) * U_im1,
        # Z_3_im1: sqrt(H_im1) * S_im1,
        # Z_4_im1: sqrt(H_im1) * K_im1,
        # Z_5_im1: sqrt(H_im1) * M_im1,
    # }


def set_generalized_shallow_water_variables_2d(num_moments):
    g = var("g", domain="positive")
    # primitive variables
    h = var("h", domain="positive")
    u = var("u")
    v = var("v")
    if num_moments > 0:
        alpha = get_vector_variable("alpha", num_moments)
        beta = get_vector_variable("beta", num_moments)

    list_ = [h, u, v]
    for i in range(num_moments):
        list_.append(alpha[i])
        list_.append(beta[i])

    global primitive, conserved
    primitive = vector(list_)
    conserved = get_vector_variable("q", 3 + 2 * num_moments)

    global p, q
    p = primitive
    q = conserved

    # substitution dictionaries
    global p_to_q, q_to_p
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
        p_to_q[alpha[i]] = q[3 + 2 * i] / q[0]
        p_to_q[beta[i]] = q[4 + 2 * i] / q[0]
        q_to_p[q[3 + 2 * i]] = h * alpha[i]
        q_to_p[q[4 + 2 * i]] = h * beta[i]

    # transformation symbolic expressions
    global p_q, q_p
    p_q = vector(p_to_q.values())
    q_p = vector(q_to_p.values())

    # transformation jacobians
    global p_q_j, q_p_j
    p_q_j = jacobian(p_q, q)
    q_p_j = jacobian(q_p, p)

    # get equations, get f_x_p, f_y_p,
    get_generalized_shallow_water_equations_2d(num_moments)

    # conserved and transform fluxes
    global f_x_q, f_y_q
    f_x_q = f_x_p.subs(p_to_q)
    f_y_q = f_y_p.subs(p_to_q)

    # flux_jacobians
    global f_x_p_j, f_x_q_j, f_y_p_j, f_y_q_j
    f_x_p_j = jacobian(f_x_p, p)
    f_y_p_j = jacobian(f_y_p, p)
    f_x_q_j = jacobian(f_x_q, q)
    f_y_q_j = jacobian(f_y_q, q)

    # Nonconserved matrix
    global G_x_q, G_y_q
    G_x_q = G_x_p.subs(p_to_q)
    G_y_q = G_y_p.subs(p_to_q)

    # Quasilinear Matrices
    global A_p, A_q, B_p, B_q, A_x_p, A_x_q, A_y_p, A_y_q
    A_x_q = f_x_q_j - G_x_q
    A_y_q = f_y_q_j - G_y_q
    # should also be equivalent to f_x_p_j @ p_q_j - G_x_p
    A_x_p = A_x_q.subs(q_to_p)
    # should also be equivalent to f_y_p_j @ p_q_j - G_y_p
    A_y_p = A_y_q.subs(q_to_p)

    A_p = A_x_p
    A_q = A_x_q
    B_p = A_y_p
    B_q = A_y_q


def get_generalized_shallow_water_equations_2d(num_moments, is_functions=False):
    # pass in number of moments not counting constant moment
    num_eqns = 3 + 2 * num_moments
    t = var("t")
    x = var("x")
    y = var("y")
    z = var("z")
    g = var("g")
    e_x = var("e_x")
    e_y = var("e_y")
    e_z = var("e_z")
    nu = var("nu")
    lambda_ = var("lambda_")
    i = var("i", domain="integer")
    j = var("j", domain="integer")
    k = var("k", domain="integer")
    h_b = function("h_b", nargs=2)(x, y)
    if is_functions:
        h = function("h", nargs=3)(t, x, y)
        u = function("u", nargs=3)(t, x, y)
        v = function("v", nargs=3)(t, x, y)
        get_symbolic = lambda str_name: function(str_name, nargs=3)(t, x, y)
        alpha = get_vector_symbolic("alpha", num_moments, get_symbolic)
        beta = get_vector_symbolic("beta", num_moments, get_symbolic)
    else:
        h = var("h", domain="positive")
        u = var("u")
        v = var("v")
        alpha = get_vector_variable("alpha", num_moments)
        beta = get_vector_variable("beta", num_moments)

    phi = get_legendre_polynomials_fixed_lower_endpoint(num_moments + 1, 1, 0, 1)
    A = [
        [
            [
                (2 * i + 3) * (phi[i + 1] * phi[j + 1] * phi[k + 1]).integrate(x, 0, 1)
                for k in range(num_moments)
            ]
            for j in range(num_moments)
        ]
        for i in range(num_moments)
    ]
    B = [
        [
            [
                (2 * i + 3)
                * (
                    phi[i + 1].derivative(x) * (phi[j + 1](x=z)).integrate(z, 0, x) * phi[k + 1]
                ).integrate(x, 0, 1)
                for k in range(num_moments)
            ]
            for j in range(num_moments)
        ]
        for i in range(num_moments)
    ]
    C = [
        [
            (phi[i + 1].derivative(x) * phi[j + 1].derivative(x)).integrate(x, 0, 1)
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
        h * u
        ^ 2
        + h * sum([alpha[j] ^ 2 / (2 * j + 3) for j in range(num_moments)])
        + 1 / 2 * g * e_z * h
        ^ 2,
        h * u * v
        + h * sum([alpha[j] * beta[j] / (2 * j + 3) for j in range(num_moments)]),
    ]
    f_y_list = [
        h * v,
        h * u * v
        + h * sum([alpha[j] * beta[j] / (2 * j + 3) for j in range(num_moments)]),
        h * v
        ^ 2
        + h * sum([beta[j] ^ 2 / (2 * j + 3) for j in range(num_moments)])
        + 1 / 2 * g * e_z * h
        ^ 2,
    ]

    # source term
    s_list = [
        0,
        -nu / lambda_ * (u + sum([alpha[j] for j in range(num_moments)]))
        + h * g * (e_x - e_z * h_b.derivative(x)),
        -nu / lambda_ * (v + sum([beta[j] for j in range(num_moments)]))
        + h * g * (e_y - e_z * h_b.derivative(y)),
    ]

    for i in range(num_moments):
        f_x_list.append(
            2 * h * u * alpha[i]
            + h
            * sum(
                [
                    sum(
                        [
                            A[i][j][k] * alpha[j] * alpha[k]
                            for k in range(num_moments)
                        ]
                    )
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
                    sum(
                        [A[i][j][k] * alpha[j] * beta[k] for k in range(num_moments)]
                    )
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
                    sum(
                        [A[i][j][k] * alpha[j] * beta[j] for k in range(num_moments)]
                    )
                    for j in range(num_moments)
                ]
            )
        )
        f_y_list.append(
            2 * h * v * beta[i]
            + h
            * sum(
                [
                    sum([A[i][j][k] * beta[j] * beta[k] for k in range(num_moments)])
                    for j in range(num_moments)
                ]
            )
        )

        s_list.append(
            -(2 * i + 3)
            * nu
            / lambda_
            * (
                u
                + sum(
                    [
                        (1 + lambda_ / h * C[i][j]) * alpha[j]
                        for j in range(num_moments)
                    ]
                )
            )
        )
        s_list.append(
            -(2 * i + 3)
            * nu
            / lambda_
            * (
                v
                + sum(
                    [
                        (1 + lambda_ / h * C[i][j]) * beta[j]
                        for j in range(num_moments)
                    ]
                )
            )
        )

    # nonconservative product matrices
    G_x_lists = [[0 for i in range(num_eqns)] for j in range(num_eqns)]
    G_y_lists = [[0 for i in range(num_eqns)] for j in range(num_eqns)]

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

    global f_x_p, f_y_p, G_x_p, G_y_p, s_p
    f_x_p = vector(f_x_list)
    f_y_p = vector(f_y_list)
    G_x_p = matrix(G_x_lists)
    G_y_p = matrix(G_y_lists)
    s_p = vector(s_list)

    return (f_x_p, f_y_p, G_x_p, G_y_p, s_p)


def get_generalized_shallow_water_equations_1d(num_moments, is_functions=False):
    if is_functions:
        v = function("v", nargs=3)(t, x, y)
        get_symbolic = lambda str_name: function(str_name, nargs=3)(t, x, y)
        beta = get_vector_symbolic("beta", num_moments, get_symbolic)
    else:
        v = var("v")
        beta = get_vector_variable("beta", num_moments)

    dict_ = {v: 0}
    for i in range(num_moments):
        dict_[beta[i]] = 0

    tuple_ = get_generalized_shallow_water_equations_2d(num_moments, is_functions)

    global f_p, G_p, s_p
    f_p = tuple_[0].subs(dict_)
    G_p = tuple_[2].subs(dict_)
    s_p = tuple_[4].subs(dict_)

    rows_delete = [2 * i + 2 for i in range(num_moments+1)]
    G_p = G_p.delete_rows(rows_delete)
    G_p = G_p.delete_columns(rows_delete)

    rows_keep = [2 * i + 1 for i in range(num_moments+1)]
    rows_keep.insert(0, 0)
    f_p = vector([f_p[i] for i in rows_keep])
    s_p = vector([s_p[i] for i in rows_keep])

    return (f_p, G_p, s_p)

