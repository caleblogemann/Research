# import all of sage underneath sage.all namespace
import sage.all as sa


def is_jacobian(A, variables):
    # check if A is the jacobian of a vector function f
    # with respect to variables
    # A is symbolic matrix
    # variables, vector of variables
    n_rows = A.nrows()
    n_cols = A.ncols()
    assert len(variables) == n_cols
    for i in range(n_rows):
        for k in range(n_cols):
            for l in range(k):
                d_a_ik_l = sa.derivative(A[i, k], variables[l])
                d_a_il_k = sa.derivative(A[i, l], variables[k])
                if ((d_a_ik_l - d_a_il_k).simplify_full() != 0 ):
                    return False

    return True


def function_from_jacobian(A, variables):
    # find function f such that jacobian of f with respect to variables is A
    # A is symbolic matrix
    # variables is symbolic vector

    num_variables = len(variables)
    assert A.ncols() == num_variables

    substitution_dict = {v : 0 for v in variables}
    f = sa.vector(sa.SR, A.nrows())
    for i_var in range(num_variables):
        substitution_dict.pop(variables[i_var])
        col = A.column(i_var).subs(substitution_dict)
        f += sa.vector([sa.integral(a, variables[i_var]) for a in col])

    return f


if __name__ == "__main__":
    pass
