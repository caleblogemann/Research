def get_vector_symbolic(name, length, get_symbolic):
    list_ = []
    for i in range(length):
        sym = get_symbolic(name + '_' + str(i))
        list_.append(sym)
    return vector(list_)


def get_matrix_symbolic(name, num_rows, num_cols, get_symbolic):
    list_1 = []
    for i in range(num_rows):
        list_2 = []
        for j in range(num_cols):
            sym = get_symbolic(name + '_' + str(i) + '_' + str(j))
            list_2.append(sym)
        list_1.append(list_2)
    return matrix(list_1)


def get_vector_variable(name, length):
    return get_vector_symbolic(name, length, var)


def get_matrix_variable(name, num_rows, num_cols):
    return get_matrix_symbolic(name, num_rows, num_cols, var)


def get_vector_single_variable_function(name, length, variable):
    get_symbolic = lambda str_name: function(str_name, nargs=1)(variable)
    return get_vector_symbolic(name, length, get_symbolic)


def get_matrix_single_variable_function(name, num_rows, num_cols, variable):
    get_symbolic = lambda str_name: function(str_name, nargs=1)(variable)
    return get_matrix_symbolic(name, num_rows, num_cols, get_symbolic)


if __name__ == "__main__":
    pass
