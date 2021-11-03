from sageresearch.utils import symbolic_vector_matrix as svm

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


class CanonicalElement:
    @staticmethod
    def get_mesh_variables():
        raise Exception()

    @staticmethod
    def get_canonical_variables():
        raise Exception()

    @staticmethod
    def get_mesh_element_vertices():
        raise Exception()

    @staticmethod
    def get_canonical_element_vertices():
        raise Exception()

    @staticmethod
    def integrate_over_mesh_element(f, vertex_list):
        raise Exception()

    @staticmethod
    def integrate_over_canonical_element(f):
        raise Exception()

    @staticmethod
    def transform_to_canonical(x, vertex_list):
        raise Exception()

    @staticmethod
    def get_transformation_to_canonical(vertex_list):
        raise Exception()

    @staticmethod
    def transform_to_canonical_jacobian(vertex_list):
        raise Exception()

    @staticmethod
    def transform_to_canonical_jacobian_determinant(vertex_list):
        raise Exception()

    @staticmethod
    def transform_to_mesh(xi, vertex_list):
        raise Exception()

    @staticmethod
    def get_transformation_to_mesh(vertex_list):
        raise Exception()

    @staticmethod
    def transform_to_mesh_jacobian(vertex_list):
        raise Exception()

    @staticmethod
    def transform_to_mesh_jacobian_determinant(vertex_list):
        raise Exception()


class Interval(CanonicalElement):
    @staticmethod
    def get_mesh_variables():
        x = sa.SR.symbol("x", domain="real")
        return x

    @staticmethod
    def get_canonical_variables():
        xi = sa.SR.symbol("xi", domain="real")
        return xi

    @staticmethod
    def get_mesh_element_vertices():
        x_l = sa.SR.symbol("x_l", domain="real")
        x_r = sa.SR.symbol("x_r", domain="real")
        return sa.matrix([[x_l], [x_r]])

    @staticmethod
    def get_canonical_element_vertices():
        return vertices_1d

    @staticmethod
    def integrate_over_mesh_element(f, vertex_list):
        # f should be symbolic expression of mesh variables, x
        x = Interval.get_mesh_variables()
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        return sa.integrate(f, x, x_left, x_right)

    @staticmethod
    def integrate_over_canonical_element(f):
        # f should be symbolic expression of canonical variables, xi
        xi = Interval.get_canonical_variables()
        return sa.integrate(f, xi, -1, 1)

    @staticmethod
    def transform_to_canonical(x, vertex_list):
        # xi = (x - x_c) 2 / delta_x
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        x_c = (x_left + x_right) / two
        delta_x = x_right - x_left
        xi = (x - x_c) * two / delta_x
        return xi

    @staticmethod
    def get_transformation_to_canonical(vertex_list):
        x = Interval.get_mesh_variables()
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        c_i = Interval.transform_to_canonical(x, vertex_list)
        return c_i.function(x)

    @staticmethod
    def transform_to_canonical_jacobian(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        delta_x = x_right - x_left
        return sa.matrix([[two / delta_x]])

    @staticmethod
    def transform_to_canonical_jacobian_determinant(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        delta_x = x_right - x_left
        return two / delta_x


    @staticmethod
    def transform_to_mesh(xi, vertex_list):
        # x = xi * delta_x / 2 + x_c
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        x_c = (x_left + x_right) / two
        delta_x = x_right - x_left
        x = xi * delta_x / two + x_c
        return x

    @staticmethod
    def get_transformation_to_mesh(vertex_list):
        # return function b_i(xi)
        xi = Interval.get_canonical_variables()
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        b_i = Interval.transform_to_mesh(xi, vertex_list)
        return b_i.function(xi)

    @staticmethod
    def transform_to_mesh_jacobian(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        delta_x = x_right - x_left
        return sa.matrix([[delta_x / two]])

    @staticmethod
    def transform_to_mesh_jacobian_determinant(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        delta_x = x_right - x_left
        return delta_x / two


class CanonicalElement2D(CanonicalElement):

    @staticmethod
    def get_mesh_variables():
        x = sa.SR.symbol("x", domain="real")
        y = sa.SR.symbol("y", domain="real")
        return (x, y)

    @staticmethod
    def get_canonical_variables():
        xi = sa.SR.symbol("xi", domain="real")
        eta = sa.SR.symbol("eta", domain="real")
        return (xi, eta)

    @staticmethod
    def get_parameterization_variable():
        s = sa.SR.symbol("s", domain="real")
        return s

    def get_transformation_to_canonical(self, vertex_list):
        tuple_ = self.get_mesh_variables()
        x = tuple_[0]
        y = tuple_[1]
        c_i = self.transform_to_canonical(tuple_, vertex_list)
        return c_i.function(x, y)

    def get_transformation_to_mesh(self, vertex_list):
        tuple_ = self.get_canonical_variables()
        xi = tuple_[0]
        eta = tuple_[1]
        b_i = self.transform_to_mesh(tuple_, vertex_list)
        return b_i.function(xi, eta)


class Square(CanonicalElement2D):
    @staticmethod
    def get_mesh_element_vertices():
        x_l = sa.SR.symbol("x_l", domain="real")
        x_r = sa.SR.symbol("x_r", domain="real")
        y_b = sa.SR.symbol("y_b", domain="real")
        y_t = sa.SR.symbol("y_t", domain="real")
        vertex_list = sa.matrix([[x_l, y_b], [x_r, y_b], [x_r, y_t], [x_l, y_t]])
        return vertex_list

    @staticmethod
    def get_canonical_element_vertices():
        return vertices_2d_rectangle

    @staticmethod
    def integrate_over_mesh_element(f, vertex_list):
        tuple_ = Square.get_mesh_variables()
        x = tuple_[0]
        y = tuple_[1]
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        integral = sa.integrate(sa.integrate(f, x, x_left, x_right), y, y_bottom, y_top)
        return integral

    @staticmethod
    def integrate_over_canonical_element(f):
        tuple_ = Square.get_canonical_variables()
        xi = tuple_[0]
        eta = tuple_[1]
        integral = sa.integrate(sa.integrate(f, xi, -one, one), eta, -one, one)
        return integral

    @staticmethod
    def transform_to_canonical(x, vertex_list):
        # vertex_list should be in order, bottom_left, bottom_right, top_right, top_left
        # xi = (x - x_c) * 2.0 / delta_x
        # eta = (y - y_c) * 2.0 / delta_x
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        x_c = (x_left + x_right) / two
        y_c = (y_bottom + y_top) / two
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        xi = (x[0] - x_c) * two / delta_x
        eta = (x[1] - y_c) * two / delta_y
        return sa.vector([xi, eta])

    @staticmethod
    def transform_to_canonical_jacobian(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        return sa.matrix([[two / delta_x, 0.0], [0, two / delta_y]])

    @staticmethod
    def transform_to_canonical_jacobian_determinant(
        vertex_list,
    ):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        return four / (delta_x * delta_y)

    @staticmethod
    def transform_to_mesh(xi, vertex_list):
        # xi should either be shape (2,) or (num_points, 2)
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        x_c = (x_left + x_right) / two
        y_c = (y_bottom + y_top) / two
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        x = xi[0] * delta_x / two + x_c
        y = xi[1] * delta_y / two + y_c
        return sa.vector([x, y])

    @staticmethod
    def transform_to_mesh_jacobian(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        return sa.matrix([[delta_x / two, 0], [0, delta_y / two]])

    @staticmethod
    def transform_to_mesh_jacobian_determinant(vertex_list):
        x_left = vertex_list[0, 0]
        x_right = vertex_list[1, 0]
        y_bottom = vertex_list[0, 1]
        y_top = vertex_list[2, 1]
        delta_x = x_right - x_left
        delta_y = y_top - y_bottom
        return (delta_x * delta_y) / four

    @staticmethod
    def get_parameterizations_of_mesh_element_faces():
        s = Square.get_parameterization_variable()
        vertex_list = Square.get_mesh_element_vertices()
        x_l = vertex_list[0, 0]
        x_r = vertex_list[1, 0]
        y_b = vertex_list[0, 1]
        y_t = vertex_list[2, 1]
        delta_y = y_t - y_b
        y_c = (y_b + y_t) / two
        delta_x = x_r - x_l
        x_c = (x_r + x_l) / two
        r_l = sa.vector([x_l, y_c + s * delta_y / two]).function(s)
        r_r = sa.vector([x_r, y_c + s * delta_y / two]).function(s)
        r_b = sa.vector([x_c + s * delta_x / two, y_b]).function(s)
        r_t = sa.vector([x_c + s * delta_x / two, y_t]).function(s)
        return (r_l, r_r, r_b, r_t)

    @staticmethod
    def get_parameterizations_of_canonical_element_faces():
        s = Square.get_parameterization_variable()
        r_l = sa.vector([-1, s]).function(s)
        r_r = sa.vector([1, s]).function(s)
        r_b = sa.vector([s, -1]).function(s)
        r_t = sa.vector([s, 1]).function(s)
        return (r_l, r_r, r_b, r_t)

    @staticmethod
    def get_outward_normal_vectors_mesh_element():
        n_l = sa.vector([-1, 0])
        n_r = sa.vector([1, 0])
        n_b = sa.vector([0, -1])
        n_t = sa.vector([0, 1])
        return (n_l, n_r, n_b, n_t)

    @staticmethod
    def get_outward_normal_vectors_canonical_element():
        return Square.get_outward_normal_vectors_mesh_element()


class Triangle(CanonicalElement2D):
    @staticmethod
    def get_mesh_element_vertices():
        x = svm.get_vector_variable("x", 3, domain="real")
        y = svm.get_vector_variable("y", 3, domain="real")
        vertex_list = sa.matrix([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]]])
        return vertex_list

    @staticmethod
    def get_canonical_element_vertices():
        return vertices_2d_triangle

    @staticmethod
    def integrate_over_mesh_element(f, vertex_list):
        tuple_ = Triangle.get_mesh_variables()
        x = tuple_[0]
        y = tuple_[1]

        sorted_vertex_list = vertex_list.copy()
        sorted_vertex_list.sort()

        x_left = sorted_vertex_list[0, 0]
        y_left = sorted_vertex_list[0, 1]
        x_middle = sorted_vertex_list[1, 0]
        y_middle = sorted_vertex_list[1, 1]
        x_right = sorted_vertex_list[2, 0]
        y_right = sorted_vertex_list[2, 1]
        slope_left_to_middle = (y_middle - y_left) / (x_middle - x_left)
        slope_left_to_right = (y_right - y_left) / (x_right - x_left)
        slope_middle_to_right = (y_right - y_middle) / (x_right - x_middle)

        y_left_to_right = slope_left_to_right * (x - x_left) + y_left
        y_left_to_middle = slope_left_to_middle * (x - x_left) + y_left
        y_middle_to_right = slope_middle_to_right * (x - x_right) + y_right

        if y_middle > y_right:
            # if middle vertex higher than right vertex
            y_bottom_left = y_left_to_right
            y_top_left = y_left_to_middle

            y_bottom_right = y_left_to_right
            y_top_right = y_middle_to_right
        else:
            # if right vertex higher than middle vertex
            y_bottom_left = y_left_to_middle
            y_top_left = y_left_to_right
            y_bottom_right = y_middle_to_right
            y_top_right = y_left_to_right

        left_integral = sa.integrate(
            sa.integrate(f, y, y_bottom_left, y_top_left), x, x_left, x_middle
        )

        right_integral = sa.integrate(
            sa.integrate(f, y, y_bottom_right, y_top_right), x, x_middle, x_right
        )

        return left_integral + right_integral

    @staticmethod
    def integrate_over_canonical_element(f):
        tuple_ = Triangle.get_canonical_variables()
        xi = tuple_[0]
        eta = tuple_[1]
        integral = sa.integrate(sa.integrate(f, xi, -one, -eta), eta, -one, one)
        return integral

    @staticmethod
    def transform_to_canonical(x, vertex_list):
        # x.shape = (2, points.shape)
        # vertex_list = (3, 2)
        # return shape (2, points.shape)
        # xi[0] = a_0_0 x[0] + a_0_1 x[1] + a_0_2
        # xi[1] = a_1_0 x[0] + a_1_1 x[1] + a_1_2
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

    @staticmethod
    def transform_to_canonical_jacobian(vertex_list):
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

    @staticmethod
    def transform_to_canonical_jacobian_determinant(vertex_list):
        # determinant of jacobian of transformation to canonical element
        # should be constant scalar as transformation is linear
        jacobian = transform_to_canonical_jacobian_2d_triangle_vertex_list(vertex_list)
        det = sa.det(jacobian)
        return det

    @staticmethod
    def transform_to_mesh(xi, vertex_list):
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

    @staticmethod
    def transform_to_mesh_jacobian(vertex_list):
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

    @staticmethod
    def transform_to_mesh_jacobian_determinant(vertex_list):
        # jacobian of transformation to mesh
        # should be constant scalar as the transformation is linear
        # x[0] = b_0_0 xi[0] + b_0_1 xi[1] + b_0_2
        # x[1] = b_1_0 xi[0] + b_1_1 xi[1] + b_1_2
        # det(J) = b_0_0 * b_1_1 - b_0_1 *b_1_0

        jacobian = transform_to_mesh_jacobian_2d_triangle_vertex_list(vertex_list)
        det = sa.det(jacobian)
        return det

    @staticmethod
    def get_parameterizations_of_mesh_element_faces():
        s = Triangle.get_parameterization_variable()
        vertex_list = Triangle.get_mesh_element_vertices()
        v_1 = vertex_list[0]
        v_2 = vertex_list[1]
        v_3 = vertex_list[2]

        # r_l v_1 to v_2
        r_l = sa.vector(((v_2 - v_1) * s + v_1 + v_2) / two).function(s)
        # r_b v_2 to v_3
        r_b = sa.vector(((v_3 - v_2) * s + v_2 + v_3) / two).function(s)
        # r_h v_3 to v_1
        r_h = sa.vector(((v_1 - v_3) * s + v_3 + v_1) / two).function(s)

        return (r_l, r_b, r_h)

    @staticmethod
    def get_parameterizations_of_canonical_element_faces():
        s = Triangle.get_parameterization_variable()
        r_l = sa.vector([-one, s]).function(s)
        r_b = sa.vector([s, -one]).function(s)
        r_h = sa.vector([s, -s]).function(s)
        return (r_l, r_b, r_h)

    @staticmethod
    def get_outward_normal_vectors_mesh_element():
        v = Triangle.get_mesh_element_vertices()

        # n_l - outward normal to face from vertex[0] to vertex[1]
        n_l = sa.vector([v[1, 1] - v[0, 1], v[0, 0] - v[1, 0]])
        n_l = n_l / sa.norm(n_l)
        # n_b - outward normal to face from vertex[1] to vertex[2]
        n_b = sa.vector([v[2, 1] - v[1, 1], v[1, 0] - v[2, 0]])
        n_b = n_b / sa.norm(n_b)
        # n_h - outward normal to face from vertex[2] to vertex[0]
        n_h = sa.vector([v[0, 1] - v[2, 1], v[2, 0] - v[0, 0]])
        n_h = n_h / sa.norm(n_h)
        return (n_l, n_b, n_h)

    @staticmethod
    def get_outward_normal_vectors_canonical_element():
        n_l = sa.vector([-one, zero])
        n_b = sa.vector([zero, -one])
        n_h = sa.vector([one / sa.sqrt(two), one / sa.sqrt(two)])
        return (n_l, n_b, n_h)
