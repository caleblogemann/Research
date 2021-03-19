import sageresearch.utils.symbolic_vector_matrix as symbolic_vector_matrix


def test_get_vector_variable():
    length = 4
    q = symbolic_vector_matrix.get_vector_variable('q', length)
    assert len(q) == length

    assert False


def test_get_matrix_variable():
    assert False


def test_get_vector_single_variable_function():
    assert False


def test_get_matrix_single_variable_function():
    assert False
