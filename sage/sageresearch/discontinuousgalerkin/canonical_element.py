import sage.all as sa

def transform_to_canonical_1d_vertex_list(x, vertex_list):
    return transform_to_canonical_1d(
        x, vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_canonical_1d(x, x_left, x_right):
    # xi = (x - x_c) 2 / delta_x
    x_c = 0.5 * (x_left + x_right)
    delta_x = x_right - x_left
    xi = (x - x_c) * 2.0 / delta_x
    return xi

def transform_to_canonical_jacobian_1d_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_canonical_jacobian_1d(x_left, x_right):
    delta_x = x_right - x_left
    return sa.matrix([[2.0 / delta_x]])

def transform_to_canonical_jacobian_determinant_1d_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_determinant_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_canonical_jacobian_determinant_1d(x_left, x_right):
    delta_x = x_right - x_left
    return 2.0 / delta_x

def transform_to_mesh_1d_vertex_list(xi, vertex_list):
    return transform_to_mesh_1d(
        xi, vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_mesh_1d(xi, x_left, x_right):
    # x = xi * delta_x / 2 + x_c
    x_c = 0.5 * (x_left + x_right)
    delta_x = x_right - x_left
    x = xi * delta_x * 0.5 + x_c
    return x

def transform_to_mesh_jacobian_1d_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_mesh_jacobian_1d(x_left, x_right):
    delta_x = x_right - x_left
    return sa.matrix([[delta_x * 0.5]])

def transform_to_mesh_jacobian_determinant_1d_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_determinant_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )

def transform_to_mesh_jacobian_determinant_1d(x_left, x_right):
    delta_x = x_right - x_left
    return 0.5 * delta_x

