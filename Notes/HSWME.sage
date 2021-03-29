def get_hswme_quasilinear_matrix(num_moments):
    num_eqns = num_moments + 2

    A = matrix(SR, (num_eqns, num_eqns))
    pass