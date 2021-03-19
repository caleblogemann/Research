import sageresearch.utils.calculus_utils as calculus_utils
import sageresearch.utils.symbolic_vector_matrix as symbolic_vector_matrix

import sage.all as sa

def test_is_jacobian():
    # test function from jacobian
    q = symbolic_vector_matrix.get_vector_variable('q', 2)
    f = sa.vector([q[0], q[1]])
    A = sa.jacobian(f, q)
    assert calculus_utils.is_jacobian(A, q)

    f = sa.vector([q[0] * q[1] + q[1], q[0] * q[1] + q[0]])
    A = sa.jacobian(f, q)
    assert calculus_utils.is_jacobian(A, q)

    f = sa.vector([sa.sin(q[1]) * sa.cos(q[0]), sa.sin(q[1] * q[0]), q[1] * sa.cos(q[0])])
    A = sa.jacobian(f, q)
    assert calculus_utils.is_jacobian(A, q)


def test_function_from_jacobian():
    # test function from jacobian
    q = symbolic_vector_matrix.get_vector_variable('q', 2)
    f = sa.vector([q[0], q[1]])
    A = sa.jacobian(f, q)
    f2 = calculus_utils.function_from_jacobian(A, f)
    assert (f - f2).norm().simplify_full() == 0.0

    f = sa.vector([q[0] * q[1] + q[1], q[0] * q[1] + q[0]])
    A = sa.jacobian(f, q)
    f2 = calculus_utils.function_from_jacobian(A, q)
    assert (f - f2).norm().simplify_full() == 0.0

    f = sa.vector([sa.sin(q[1]) * sa.cos(q[0]), sa.sin(q[1] * q[0]), q[1] * sa.cos(q[0])])
    A = sa.jacobian(f, q)
    f2 = calculus_utils.function_from_jacobian(A, q)
    assert (f - f2).norm().simplify_full() == 0.0

