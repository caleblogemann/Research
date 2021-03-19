import sageresearch.shallowwater.shallow_water_moment_equations as swme

import sage.all as sa

max_num_moments = 3

def test_get_conserved_variables_1d():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        q = swme.get_conserved_variables_1d(num_moments)
        assert len(q) == num_eqns


def test_get_primitive_variables_1d():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        p = swme.get_primitive_variables_1d(num_moments)
        assert len(p) == num_eqns


def test_get_substitution_dictionaries_1d():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        tuple_ = swme.get_substitution_dictionaries_1d(num_moments)
        p_to_q = tuple_[0]
        q_to_p = tuple_[1]
        assert len(p_to_q) == num_eqns
        assert len(q_to_p) == num_eqns
        p = swme.get_primitive_variables_1d(num_moments)
        assert sa.vector(p_to_q.values()).subs(q_to_p) == p
        q = swme.get_conserved_variables_1d(num_moments)
        assert sa.vector(q_to_p.values()).subs(p_to_q) == q


def test_get_conserved_variables_2d():
    for num_moments in range(max_num_moments):
        num_eqns = 2 * num_moments + 3
        q = swme.get_conserved_variables_2d(num_moments)
        assert len(q) == num_eqns


def test_get_primitive_variables_2d():
    for num_moments in range(max_num_moments):
        num_eqns = 2 * num_moments + 3
        p = swme.get_primitive_variables_2d(num_moments)
        assert len(p) == num_eqns


def test_get_substitution_dictionaries_2d():
    for num_moments in range(max_num_moments):
        num_eqns = 2 * num_moments + 3
        tuple_ = swme.get_substitution_dictionaries_2d(num_moments)
        p_to_q = tuple_[0]
        q_to_p = tuple_[1]
        assert len(p_to_q) == num_eqns
        assert len(q_to_p) == num_eqns
        p = swme.get_primitive_variables_2d(num_moments)
        assert sa.vector(p_to_q.values()).subs(q_to_p) == p
        q = swme.get_conserved_variables_2d(num_moments)
        assert sa.vector(q_to_p.values()).subs(p_to_q) == q


def test_get_generalized_shallow_water_equations_1d():
    for num_moments in range(max_num_moments):
        num_eqns = num_moments + 2
        tuple_ = swme.get_generalized_shallow_water_equations_1d(num_moments)
        f_p = tuple_[0]
        G_p = tuple_[1]
        s_p = tuple_[2]
        A_p = tuple_[3]

        assert len(f_p) == num_eqns
        assert G_p.nrows() == num_eqns
        assert G_p.ncols() == num_eqns
        assert len(s_p) == num_eqns
        assert A_p.nrows() == num_eqns
        assert A_p.ncols() == num_eqns


def test_get_generalized_shallow_water_equations_2d():
    for num_moments in range(max_num_moments):
        num_eqns = 2 * num_moments + 3
        tuple_ = swme.get_generalized_shallow_water_equations_1d(num_moments)
        f_x_p = tuple_[0]
        f_y_p = tuple_[1]
        G_x_p = tuple_[2]
        G_y_p = tuple_[3]
        s_p = tuple_[4]
        A_x_p = tuple_[5]
        A_y_p = tuple_[6]

        assert len(f_x_p) == num_eqns
        assert len(f_y_p) == num_eqns
        assert G_x_p.nrows() == num_eqns
        assert G_x_p.ncols() == num_eqns
        assert G_y_p.nrows() == num_eqns
        assert G_y_p.ncols() == num_eqns
        assert len(s_p) == num_eqns
        assert A_x_p.nrows() == num_eqns
        assert A_x_p.ncols() == num_eqns
        assert A_y_p.nrows() == num_eqns
        assert A_y_p.ncols() == num_eqns
