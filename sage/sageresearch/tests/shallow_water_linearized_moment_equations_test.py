import sageresearch.shallowwater.shallow_water_linearized_moment_equations as swlme
import sageresearch.shallowwater.shallow_water_moment_equations as swme


def test_get_swlme_equations():
    for num_moments in range(5):
        tuple_ = swlme.get_swlme_equations(num_moments)
        f = tuple_[0]
        G = tuple_[1]
        A = tuple_[2]

        num_eqns = num_moments + 2
        assert len(f) == num_eqns
        assert G.nrows() == num_eqns
        assert G.ncols() == num_eqns
        assert A.nrows() == num_eqns
        assert A.ncols() == num_eqns


def test_swlme_equal_to_swme():
    # for 0 and 1 moments swlme should be equivalent to standard shallow water equations
    for num_moments in range(2):
        tuple_ = swlme.get_swlme_equations(num_moments)
        f_swlme = tuple_[0]
        G_swlme = tuple_[1]
        A_swlme = tuple_[2]

        tuple_ = swme.get_generalized_shallow_water_equations_1d(num_moments)
        f_swme = tuple_[0]
        G_swme = tuple_[1]
        A_swme = tuple_[2]

        assert f_swlme == f_swme
        assert G_swlme == G_swme
        assert A_swlme == A_swme

