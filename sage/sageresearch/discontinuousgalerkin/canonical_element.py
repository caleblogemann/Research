import sage.all as sa

zero = sa.Integer(0)
one = sa.Integer(1)
two = sa.Integer(2)
four = sa.Integer(4)

volume_1d = two
volume_2d_rectangle = four
volume_2d_triangle = two

vertices_1d = sa.matrix([[-one], [one]])
vertices_2d_rectangle = sa.matrix([[-one, -one], [one, -one], [one, one], [-one, one]])
vertices_2d_triangle = sa.matrix([[-one, one], [-one, -one], [one, -one]])


def get_mesh_variables_1d():
    x = sa.SR.symbol("x", domain="real")
    return x


def get_canonical_variables_1d():
    xi = sa.SR.symbol("xi", domain="real")
    return xi


def integrate_over_mesh_element_1d(f, x_left, x_right):
    # f should be symbolic expression of mesh variables, x
    x = get_mesh_variables_1d()[0]
    return sa.integrate(f, x, x_left, x_right)


def integrate_over_canonical_element_1d(f, var):
    # f should be symbolic expression of canonical variables, xi
    xi = get_canonical_variables_1d()
    return sa.integrate(f, xi, -1, 1)


def transform_to_canonical_1d_vertex_list(x, vertex_list):
    return transform_to_canonical_1d(x, vertex_list[0, 0], vertex_list[1, 0])


def transform_to_canonical_1d(x, x_left, x_right):
    # xi = (x - x_c) 2 / delta_x
    x_c = (x_left + x_right) / two
    delta_x = x_right - x_left
    xi = (x - x_c) * two / delta_x
    return xi


def transform_to_canonical_jacobian_1d_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_1d(vertex_list[0, 0], vertex_list[1, 0])


def transform_to_canonical_jacobian_1d(x_left, x_right):
    delta_x = x_right - x_left
    return sa.matrix([[two / delta_x]])


def transform_to_canonical_jacobian_determinant_1d_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_determinant_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )


def transform_to_canonical_jacobian_determinant_1d(x_left, x_right):
    delta_x = x_right - x_left
    return two / delta_x


def transform_to_mesh_1d_vertex_list(xi, vertex_list):
    return transform_to_mesh_1d(xi, vertex_list[0, 0], vertex_list[1, 0])


def transform_to_mesh_1d(xi, x_left, x_right):
    # x = xi * delta_x / 2 + x_c
    x_c = (x_left + x_right) / two
    delta_x = x_right - x_left
    x = xi * delta_x / two + x_c
    return x


def transform_to_mesh_jacobian_1d_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_1d(vertex_list[0, 0], vertex_list[1, 0])


def transform_to_mesh_jacobian_1d(x_left, x_right):
    delta_x = x_right - x_left
    return sa.matrix([[delta_x / two]])


def transform_to_mesh_jacobian_determinant_1d_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_determinant_1d(
        vertex_list[0, 0], vertex_list[1, 0]
    )


def transform_to_mesh_jacobian_determinant_1d(x_left, x_right):
    delta_x = x_right - x_left
    return delta_x / two


def get_mesh_variables_2d():
    x = sa.SR.symbol("x", domain="real")
    y = sa.SR.symbol("y", domain="real")
    return (x, y)


def get_canonical_variables_2d():
    xi = sa.SR.symbol("xi", domain="real")
    eta = sa.SR.symbol("eta", domain="real")
    return (xi, eta)


def integrate_over_mesh_element_2d_rectangle(f, x_left, x_right, y_bottom, y_top):
    tuple_ = get_mesh_variables_2d()
    x = tuple_[0]
    y = tuple_[1]
    integral = sa.integrate(sa.integrate(f, x, x_left, x_right), y, y_bottom, y_top)
    return integral

def integrate_over_canonical_element_2d_rectangle(f):
    tuple_ = get_canonical_variables_2d()
    xi = tuple_[0]
    eta = tuple_[1]
    integral = sa.integrate(sa.integrate(f, xi, -one, one), eta, -one, one)
    return integral

def transform_to_canonical_2d_rectangle_vertex_list(x, vertex_list):
    # vertex_list should be in order, bottom_left, bottom_right, top_right, top_left
    return transform_to_canonical_2d_rectangle_interval(
        x, vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1],
    )


def transform_to_canonical_2d_rectangle_interval(x, x_left, x_right, y_bottom, y_top):
    # x should either be shape (2,) or (2, num_points)
    # xi = (x - x_c) * 2.0 / delta_x
    # eta = (y - y_c) * 2.0 / delta_x
    x_c = (x_left + x_right) / two
    y_c = (y_bottom + y_top) / two
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    xi = (x[0] - x_c) * two / delta_x
    eta = (x[1] - y_c) * two / delta_y
    return sa.vector([xi, eta])


def transform_to_canonical_jacobian_2d_rectangle_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_2d_rectangle_interval(
        vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1]
    )


def transform_to_canonical_jacobian_2d_rectangle_interval(
    x_left, x_right, y_bottom, y_top
):
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    return sa.matrix([[two / delta_x, 0.0], [0, two / delta_y]])


def transform_to_canonical_jacobian_determinant_2d_rectangle_vertex_list(vertex_list):
    return transform_to_canonical_jacobian_determinant_2d_rectangle_interval(
        vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1]
    )


def transform_to_canonical_jacobian_determinant_2d_rectangle_interval(
    x_left, x_right, y_bottom, y_top
):
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    return four / (delta_x * delta_y)


def transform_to_mesh_2d_rectangle_vertex_list(xi, vertex_list):
    return transform_to_mesh_2d_rectangle_interval(
        xi, vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1],
    )


def transform_to_mesh_2d_rectangle_interval(xi, x_left, x_right, y_bottom, y_top):
    # xi should either be shape (2,) or (num_points, 2)
    x_c = (x_left + x_right) / two
    y_c = (y_bottom + y_top) / two
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    x = xi[0] * delta_x / two + x_c
    y = xi[1] * delta_y / two + y_c
    return sa.vector([x, y])


def transform_to_mesh_jacobian_2d_rectangle_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_2d_rectangle_interval(
        vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1]
    )


def transform_to_mesh_jacobian_2d_rectangle_interval(x_left, x_right, y_bottom, y_top):
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    return sa.matrix([[delta_x / two, 0], [0, delta_y / two]])


def transform_to_mesh_jacobian_determinant_2d_rectangle_vertex_list(vertex_list):
    return transform_to_mesh_jacobian_determinant_2d_rectangle_interval(
        vertex_list[0, 0], vertex_list[1, 0], vertex_list[0, 1], vertex_list[3, 1]
    )


def transform_to_mesh_jacobian_determinant_2d_rectangle_interval(
    x_left, x_right, y_bottom, y_top
):
    delta_x = x_right - x_left
    delta_y = y_top - y_bottom
    return (delta_x * delta_y) / four


def integrate_over_mesh_element_2d_triangle(f, vertex_list):
    tuple_ = get_mesh_variables_2d()
    x = tuple_[0]
    y = tuple_[1]

    sorted_vertex_list = vertex_list.copy()
    sorted_vertex_list.sort()

    x_left = sorted_vertex_list[0, 0]
    x_middle = sorted_vertex_list[1, 0]
    x_right = sorted_vertex_list[2, 0]
    # if middle vertex higher than right vertex
    if sorted_vertex_list[1, 1] > sorted_vertex_list[2, 1]:
        y_bottom =
        y_top =
    else:

    left_integral = sa.integrate(sa.integrate(f, y, ), x, x_left, x_middle)
    right_integral = sa.integrate(sa.integrate(f, y, ), x, x_middle, x_right)

    return left_integral + right_integral



def integrate_over_canonical_element_2d_triangle(f):
    tuple_ = get_canonical_variables_2d()
    xi = tuple_[0]
    eta = tuple_[1]
    integral = sa.integrate(sa.integrate(f, xi, -one, -eta), eta, -one, one)
    return integral


def transform_to_canonical_2d_triangle_vertex_list(x, vertex_list):
    # x.shape = (2, points.shape)
    # vertex_list = (3, 2)
    # return shape (2, points.shape)
    # xi[0] = a_0_0 x[0] + a_0_1 x[1] + a_0_2
    # xi[1] = a_1_0 x[0] + a_1_1 x[1] + a_1_2
    x_0 = vertex_list[0]
    x_1 = vertex_list[1]
    x_2 = vertex_list[2]
    a_0_0 = (
        -2.0
        * (x_0[1] - x_1[1])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_0_1 = (
        2.0
        * (x_0[0] - x_1[0])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_0_2 = (
        x_0[1] * (x_1[0] + x_2[0])
        - x_0[0] * (x_1[1] + x_2[1])
        - x_1[1] * x_2[0]
        + x_1[0] * x_2[1]
    ) / (
        x_0[1] * (x_1[0] - x_2[0])
        - x_0[0] * (x_1[1] - x_2[1])
        + x_1[1] * x_2[0]
        - x_1[0] * x_2[1]
    )
    a_1_0 = (
        -2.0
        * (x_1[1] - x_2[1])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_1_1 = (
        2.0
        * (x_1[0] - x_2[0])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_1_2 = -(
        x_0[1] * (x_1[0] - x_2[0])
        - x_0[0] * (x_1[1] - x_2[1])
        - x_1[1] * x_2[0]
        + x_1[0] * x_2[1]
    ) / (
        x_0[1] * (x_1[0] - x_2[0])
        - x_0[0] * (x_1[1] - x_2[1])
        + x_1[1] * x_2[0]
        - x_1[0] * x_2[1]
    )

    xi = sa.vector(
        [a_0_0 * x[0] + a_0_1 * x[1] + a_0_2, a_1_0 * x[0] + a_1_1 * x[1] + a_1_2]
    )
    return xi


def transform_to_canonical_jacobian_2d_triangle_vertex_list(vertex_list):
    # jacobian of transformation to canonical
    # should be a constant matrix as transformation is linear
    # xi[0] = a_0_0 x[0] + a_0_1 x[1] + a_0_2
    # xi[1] = a_1_0 x[0] + a_1_1 x[1] + a_1_2
    # J = [[a_0_0, a_0_1], [a_1_0, a_1_1]]
    x_0 = vertex_list[0]
    x_1 = vertex_list[1]
    x_2 = vertex_list[2]
    a_0_0 = (
        -two
        * (x_0[1] - x_1[1])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_0_1 = (
        two
        * (x_0[0] - x_1[0])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_1_0 = (
        -two
        * (x_1[1] - x_2[1])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )
    a_1_1 = (
        two
        * (x_1[0] - x_2[0])
        / (
            x_0[1] * (x_1[0] - x_2[0])
            - x_0[0] * (x_1[1] - x_2[1])
            + x_1[1] * x_2[0]
            - x_1[0] * x_2[1]
        )
    )

    jacobian = sa.matrix([[a_0_0, a_0_1], [a_1_0, a_1_1]])
    return jacobian


def transform_to_canonical_jacobian_determinant_2d_triangle_vertex_list(vertex_list):
    # determinant of jacobian of transformation to canonical element
    # should be constant scalar as transformation is linear
    jacobian = transform_to_canonical_jacobian_2d_triangle_vertex_list(vertex_list)
    det = sa.det(jacobian)
    return det


def transform_to_mesh_2d_triangle_vertex_list(xi, vertex_list):
    # xi.shape (2, points.shape)
    # vertex_list.shape = (3, 2)
    # return shape (2, points.shape)

    # x[0] = b_0_0 xi[0] + b_0_1 xi[1] + b_0_2
    # x[1] = b_1_0 xi[0] + b_1_1 xi[1] + b_1_2
    x_0 = vertex_list[0]
    x_1 = vertex_list[1]
    x_2 = vertex_list[2]
    b_0_0 = x_2[0] / two - x_1[0] / two
    b_0_1 = x_0[0] / two - x_1[0] / two
    b_0_2 = x_0[0] / two + x_2[0] / two
    b_1_0 = x_2[1] / two - x_1[1] / two
    b_1_1 = x_0[1] / two - x_1[1] / two
    b_1_2 = x_0[1] / two + x_2[1] / two
    x = sa.vector(
        [b_0_0 * xi[0] + b_0_1 * xi[1] + b_0_2, b_1_0 * xi[0] + b_1_1 * xi[1] + b_1_2,]
    )
    return x


def transform_to_mesh_jacobian_2d_triangle_vertex_list(vertex_list):
    # jacobian of transformation to mesh
    # should be constant matrix as the transformation is linear
    # x[0] = b_0_0 xi[0] + b_0_1 xi[1] + b_0_2
    # x[1] = b_1_0 xi[0] + b_1_1 xi[1] + b_1_2
    # J = [[b_0_0, b_0_1], [b_1_0, b_1_1]]
    x_0 = vertex_list[0]
    x_1 = vertex_list[1]
    x_2 = vertex_list[2]
    b_0_0 = x_2[0] / two - x_1[0] / two
    b_0_1 = x_0[0] / two - x_1[0] / two
    b_1_0 = x_2[1] / two - x_1[1] / two
    b_1_1 = x_0[1] / two - x_1[1] / two
    jacobian = sa.matrix([[b_0_0, b_0_1], [b_1_0, b_1_1]])
    return jacobian


def transform_to_mesh_jacobian_determinant_2d_triangle_vertex_list(vertex_list):
    # jacobian of transformation to mesh
    # should be constant scalar as the transformation is linear
    # x[0] = b_0_0 xi[0] + b_0_1 xi[1] + b_0_2
    # x[1] = b_1_0 xi[0] + b_1_1 xi[1] + b_1_2
    # det(J) = b_0_0 * b_1_1 - b_0_1 *b_1_0

    jacobian = transform_to_mesh_jacobian_2d_triangle_vertex_list(vertex_list)
    det = sa.det(jacobian)
    return det
