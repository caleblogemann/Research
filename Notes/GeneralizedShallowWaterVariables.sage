def set_generalized_shallow_water_variables(num_moments):
    g = var("g")
    # primitive variables
    h = var("h")
    u = var("u")
    s = var("s")
    k = var("k")
    m = var("m")
    primitive_0 = vector([h, u])
    primitive_1 = vector([h, u, s])
    primitive_2 = vector([h, u, s, k])
    primitive_3 = vector([h, u, s, k, m])

    # f_p(p)
    flux_primitive_0 = vector([h * u, h * u ^ 2 + 1 / 2 * g * h ^ 2])
    flux_primitive_1 = vector(
        [h * u, h * u ^ 2 + 1 / 2 * g * h ^ 2 + 1 / 3 * h * s ^ 2, 2 * h * u * s]
    )
    flux_primitive_2 = vector(
        [
            h * u,
            h * u ^ 2 + 1 / 2 * g * h ^ 2 + 1 / 3 * h * s ^ 2 + 1 / 5 * h * k ^ 2,
            2 * h * u * s + 4 / 5 * h * s * k,
            2 * h * u * k + 2 / 3 * h * s ^ 2 + 2 / 7 * h * k ^ 2,
        ]
    )
    flux_primitive_3 = vector(
        [
            h * u,
            h * u
            ^ 2 + 1 / 2 * g * h
            ^ 2 + 1 / 3 * h * s
            ^ 2 + 1 / 5 * h * k
            ^ 2 + 1 / 7 * h * m
            ^ 2,
            2 * h * u * s + 4 / 5 * h * s * k + 18 / 35 * h * k * m,
            2 * h * u * k + 2 / 3 * h * s
            ^ 2 + 2 / 7 * h * k
            ^ 2 + 4 / 21 * h * m
            ^ 2 + 6 / 7 * h * s * m,
            2 * h * u * m + 6 / 5 * h * s * k + 8 / 15 * h * k * m,
        ]
    )

    # conserved variables
    q_1 = var("q_1")
    q_2 = var("q_2")
    q_3 = var("q_3")
    q_4 = var("q_4")
    q_5 = var("q_5")
    # q = [q_1, q_2]
    conserved_0 = vector([q_1, q_2])
    conserved_1 = vector([q_1, q_2, q_3])
    conserved_2 = vector([q_1, q_2, q_3, q_4])
    conserved_3 = vector([q_1, q_2, q_3, q_4, q_5])

    # transform variables
    z_1 = var("z_1")
    z_2 = var("z_2")
    z_3 = var("z_3")
    z_4 = var("z_4")
    z_5 = var("z_5")
    transform_0 = vector([z_1, z_2])
    transform_1 = vector([z_1, z_2, z_3])
    transform_2 = vector([z_1, z_2, z_3, z_4])
    transform_3 = vector([z_1, z_2, z_3, z_4, z_5])

    global primitive, p, conserved, q, transform, z, flux_primitive, f_p
    if num_moments == 0:
        primitive = primitive_0
        p = primitive
        conserved = conserved_0
        q = conserved
        transform = transform_0
        z = transform
        flux_primitive = flux_primitive_0
        f_p = flux_primitive
    elif num_moments == 1:
        primitive = primitive_1
        p = primitive
        conserved = conserved_1
        q = conserved
        transform = transform_1
        z = transform
        flux_primitive = flux_primitive_1
        f_p = flux_primitive
    elif num_moments == 2:
        primitive = primitive_2
        p = primitive
        conserved = conserved_2
        q = conserved
        transform = transform_2
        z = transform
        flux_primitive = flux_primitive_2
        f_p = flux_primitive
    elif num_moments == 3:
        primitive = primitive_3
        p = primitive
        conserved = conserved_3
        q = conserved
        transform = transform_3
        z = transform
        flux_primitive = flux_primitive_3
        f_p = flux_primitive

    # substitution dictionaries
    global p_to_q, p_to_z, q_to_p, q_to_z, z_to_p, z_to_q
    primitive_to_conserved = {
        h: q_1,
        u: q_2 / q_1,
        s: q_3 / q_1,
        k: q_4 / q_1,
        m: q_5 / q_1,
    }
    p_to_q = primitive_to_conserved
    primitive_to_transform = {
        h: z_1 ^ 2,
        u: z_2 / z_1,
        s: z_3 / z_1,
        k: z_4 / z_1,
        m: z_5 / z_1,
    }
    p_to_z = primitive_to_transform
    conserved_to_primitive = {q_1: h, q_2: h * u, q_3: h * s, q_4: h * k, q_5: h * m}
    q_to_p = conserved_to_primitive
    conserved_to_transform = {
        q_1: z_1 ^ 2,
        q_2: z_1 * z_2,
        q_3: z_1 * z_3,
        q_4: z_1 * z_4,
        q_5: z_1 * z_5,
    }
    q_to_z = conserved_to_transform
    transform_to_primitive = {
        z_1: sqrt(h),
        z_2: sqrt(h) * u,
        z_3: sqrt(h) * s,
        z_4: sqrt(h) * k,
        z_5: sqrt(h) * m,
    }
    z_to_p = transform_to_primitive
    transform_to_conserved = {
        z_1: sqrt(q_1),
        z_2: q_2 / sqrt(q_1),
        z_3: q_3 / sqrt(q_1),
        z_4: q_4 / sqrt(q_1),
        z_5: q_5 / sqrt(q_1),
    }
    z_to_q = transform_to_conserved

    # transformation symbolic expressions
    global p_q, p_z, q_p, q_z, z_p, z_q
    # p(q)
    temp = vector([q_1, q_2 / q_1, q_3 / q_1, q_4 / q_1, q_5 / q_1])
    primitive_from_conserved = temp[0 : (num_moments + 2)]
    p_q = primitive_from_conserved
    # p(z)
    temp = vector([z_1 ^ 2, z_2 / z_1, z_3 / z_1, z_4 / z_1, z_5 / z_1])
    primitive_from_transform = temp[0 : (num_moments + 2)]
    p_z = primitive_from_transform

    # q(p)
    temp = vector([h, h * u, h * s, h * k, h * m])
    conserved_from_primitive = temp[0 : (num_moments + 2)]
    q_p = conserved_from_primitive
    # q(z)
    temp = vector([z_1 ^ 2, z_1 * z_2, z_1 * z_3, z_1 * z_4, z_1 * z_4])
    conserved_from_transform = temp[0 : (num_moments + 2)]
    q_z = conserved_from_transform

    # z(p)
    temp = vector([sqrt(h), sqrt(h) * u, sqrt(h) * s, sqrt(h) * k, sqrt(h) * m])
    transform_from_primitive = temp[0 : (num_moments + 2)]
    z_p = transform_from_primitive
    # z(q)
    temp = vector(
        [sqrt(q_1), q_2 / sqrt(q_1), q_3 / sqrt(q_1), q_4 / sqrt(q_1), q_5 / sqrt(q_1)]
    )
    transform_from_conserved = temp[0 : (num_moments + 2)]
    z_q = transform_from_conserved

    # transformation jacobians
    global p_q_j, p_z_j, q_p_j, q_z_j, z_p_j, z_q_j
    p_q_j = jacobian(p_q, q)
    p_z_j = jacobian(p_z, z)
    q_p_j = jacobian(q_p, p)
    q_z_j = jacobian(q_z, z)
    z_p_j = jacobian(z_p, p)
    z_q_j = jacobian(z_q, q)

    # conserved and transform fluxes
    global flux_conserved, f_q, flux_transform, f_z
    flux_conserved = flux_primitive.subs(p_to_q)
    f_q = flux_conserved
    flux_transform = flux_primitive.subs(p_to_z)
    f_z = flux_transform

    # flux_jacobians
    global f_p_j, f_q_j, f_z_j
    global flux_primitive_jacobian, flux_conserved_jacobian, flux_transform_jacobian
    flux_primitive_jacobian = jacobian(f_p, p)
    f_p_j = flux_primitive_jacobian
    flux_conserved_jacobian = jacobian(f_q, q)
    f_q_j = flux_conserved_jacobian
    flux_transform_jacobian = jacobian(f_z, z)
    f_z_j = flux_transform_jacobian

    # Riemann Problem
    global p_to_P_i, P_i, p_to_P_im1, P_im1
    # P_i = [H_i, U_i]
    H_i = var("H_i")
    U_i = var("U_i")
    S_i = var("S_i")
    K_i = var("K_i")
    M_i = var("M_i")
    p_to_P_i = {h: H_i, u: U_i, s: S_i, k: K_i, m: M_i}
    P_i = p.subs(p_to_P_i)
    # P_im1 = [H_im1, U_im1]
    H_im1 = var("H_im1")
    U_im1 = var("U_im1")
    S_im1 = var("S_im1")
    K_im1 = var("K_im1")
    M_im1 = var("M_im1")
    p_to_P_im1 = {h: H_im1, u: U_im1, s: S_im1}
    P_im1 = p.subs(p_to_P_im1)
    global q_to_Q_i, Q_i, q_to_Q_im1, Q_im1
    # Q_i = [Q_1_i, Q_2_i]
    Q_1_i = var("Q_1_i")
    Q_2_i = var("Q_2_i")
    Q_3_i = var("Q_3_i")
    Q_4_i = var("Q_4_i")
    Q_5_i = var("Q_5_i")
    q_to_Q_i = {q_1: Q_1_i, q_2: Q_2_i, q_3: Q_3_i, q_4: Q_4_i, q_5: Q_5_i}
    Q_i = q.subs(q_to_Q_i)
    # Q_im1 = [Q_1_im1, Q_2_im1]
    Q_1_im1 = var("Q_1_im1")
    Q_2_im1 = var("Q_2_im1")
    Q_3_im1 = var("Q_3_im1")
    Q_4_im1 = var("Q_4_im1")
    Q_5_im1 = var("Q_5_im1")
    q_to_Q_im1 = {q_1: Q_1_im1, q_2: Q_2_im1, q_3: Q_3_im1, q_4: Q_4_im1, q_5: Q_5_im1}
    Q_im1 = q.subs(q_to_Q_im1)
    global z_to_Z_i, Z_i, z_to_Z_im1, Z_im1
    # Z_i = [Z_1_i, Z_2_i]
    Z_1_i = var("Z_1_i")
    Z_2_i = var("Z_2_i")
    Z_3_i = var("Z_3_i")
    Z_4_i = var("Z_4_i")
    Z_5_i = var("Z_5_i")
    z_to_Z_i = {z_1: Z_1_i, z_2: Z_2_i, z_3: Z_3_i, z_4: Z_4_i, z_5: Z_5_i}
    Z_i = z.subs(z_to_Z_i)
    # Z_im1 = [Z_1_im1, Z_2_im1]
    Z_1_im1 = var("Z_1_im1")
    Z_2_im1 = var("Z_2_im1")
    Z_3_im1 = var("Z_3_im1")
    Z_4_im1 = var("Z_4_im1")
    Z_5_im1 = var("Z_5_im1")
    z_to_Z_im1 = {z_1: Z_1_im1, z_2: Z_2_im1, z_3: Z_3_im1, z_4: Z_4_im1, z_5: Z_5_im1}
    Z_im1 = z.subs(z_to_Z_im1)

    global P_i_to_Q_i, P_i_to_Z_i
    global Q_i_to_P_i, Q_i_to_Z_i
    global Z_i_to_P_i, Z_i_to_Q_i
    P_i_to_Q_i = {
        H_i: Q_1_i,
        U_i: Q_2_i / Q_1_i,
        S_i: Q_3_i / Q_1_i,
        K_i: Q_4_i / Q_1_i,
        M_i: Q_5_i / Q_1_i,
    }
    P_i_to_Z_i = {
        H_i: Z_1_i ^ 2,
        U_i: Z_2_i / Z_1_i,
        S_i: Z_3_i / Z_1_i,
        K_i: Z_4_i / Z_1_i,
        M_i: Z_5_i / Z_1_i,
    }
    Q_i_to_P_i = {
        Q_1_i: H_i,
        Q_2_i: H_i * U_i,
        Q_3_i: H_i * S_i,
        Q_4_i: H_i * K_i,
        Q_5_i: H_i * M_i,
    }
    Q_i_to_Z_i = {
        Q_1_i: Z_1_i ^ 2,
        Q_2_i: Z_1_i * Z_2_i,
        Q_3_i: Z_1_i * Z_3_i,
        Q_4_i: Z_1_i * Z_4_i,
        Q_5_i: Z_1_i * Z_5_i,
    }
    Z_i_to_Q_i = {
        Z_1_i: sqrt(Q_1_i),
        Z_2_i: Q_2_i / sqrt(Q_1_i),
        Z_3_i: Q_3_i / sqrt(Q_1_i),
        Z_4_i: Q_4_i / sqrt(Q_1_i),
        Z_5_i: Q_5_i / sqrt(Q_1_i),
    }
    Z_i_to_P_i = {
        Z_1_i: sqrt(H_i),
        Z_2_i: sqrt(H_i) * U_i,
        Z_3_i: sqrt(H_i) * S_i,
        Z_4_i: sqrt(H_i) * K_i,
        Z_5_i: sqrt(H_i) * M_i,
    }

    global P_im1_to_Q_im1, P_im1_to_Z_im1
    global Q_im1_to_P_im1, Q_im1_to_Z_im1
    global Z_im1_to_P_im1, Z_im1_to_Q_im1
    P_im1_to_Q_im1 = {
        H_im1: Q_1_im1,
        U_im1: Q_2_im1 / Q_1_im1,
        S_im1: Q_3_im1 / Q_1_im1,
        K_im1: Q_4_im1 / Q_1_im1,
        M_im1: Q_5_im1 / Q_1_im1,
    }
    P_im1_to_Z_im1 = {
        H_im1: Z_1_im1 ^ 2,
        U_im1: Z_2_im1 / Z_1_im1,
        S_im1: Z_3_im1 / Z_1_im1,
        K_im1: Z_4_im1 / Z_1_im1,
        M_im1: Z_5_im1 / Z_1_im1,
    }
    Q_im1_to_P_im1 = {
        Q_1_im1: H_im1,
        Q_2_im1: H_im1 * U_im1,
        Q_3_im1: H_im1 * S_im1,
        Q_4_im1: H_im1 * K_im1,
        Q_5_im1: H_im1 * M_im1,
    }
    Q_im1_to_Z_im1 = {
        Q_1_im1: Z_1_im1 ^ 2,
        Q_2_im1: Z_1_im1 * Z_2_im1,
        Q_3_im1: Z_1_im1 * Z_3_im1,
        Q_4_im1: Z_1_im1 * Z_4_im1,
        Q_5_im1: Z_1_im1 * Z_5_im1,
    }
    Z_im1_to_Q_im1 = {
        Z_1_im1: sqrt(Q_1_im1),
        Z_2_im1: Q_2_im1 / sqrt(Q_1_im1),
        Z_3_im1: Q_3_im1 / sqrt(Q_1_im1),
        Z_4_im1: Q_4_im1 / sqrt(Q_1_im1),
        Z_5_im1: Q_5_im1 / sqrt(Q_1_im1),
    }
    Z_im1_to_P_im1 = {
        Z_1_im1: sqrt(H_im1),
        Z_2_im1: sqrt(H_im1) * U_im1,
        Z_3_im1: sqrt(H_im1) * S_im1,
        Z_4_im1: sqrt(H_im1) * K_im1,
        Z_5_im1: sqrt(H_im1) * M_im1,
    }

