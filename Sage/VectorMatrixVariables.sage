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
