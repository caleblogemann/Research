def get_vector_variable(name, length):
    list_ = []
    for i in range(length):
        variable = var(name + '_' + str(i))
        list_.append(variable)
    return vector(list_)


def get_matrix_variable(name, num_rows, num_cols):
    list_ = []
    for i in range(num_rows):
        list_2 = []
        for j in range(num_cols):
            variable = var(name + '_' + str(i) + '_' + str(j))
            list_2.append(variable)
        list_.append(list_2)
    return matrix(list_)


def get_vector_function(name, length, nargs):
    list_ = []
    for i in range(length):
        func = function(name + '_' + str(i), nargs=nargs)
        list_.append(func)
    return vector(list_)


def get_matrix_function(name, num_rows, num_cols, nargs):
    list_1 = []
    for i in range(num_rows):
        list_2 = []
        for j in range(num_cols):
            func = function(name + '_' + str(i) + '_' + str(j), nargs=nargs)
            list_2.append(func)
        list_1.append(list_2)
    return matrix(list_1)
