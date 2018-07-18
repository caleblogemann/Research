function [Aq] = ThinFilmDiffusionDiagonalOperator(q, dx, original_q, bc)
    [num_elems, num_basis_cpts] = size(q);

    % number of elements between cells that don't interact
    % i.e. a cell is influenced by the 2 cells to the left and the 2 cells
    % to the right of the current cell
    num_elem_spread = 3;

    Aq = zeros(num_elems, num_basis_cpts);
    qTemp = zeros(num_elems, num_basis_cpts);

    for i = 1:num_elem_spread
        for k = 1:num_basis_cpts
            qTemp(i:num_elem_spread:num_elems, k) = q(i:num_elem_spread:num_elems, k);

            AqTemp = ThinFilmDiffusionOperator(qTemp, dx, original_q, bc);
            Aq(i:num_elem_spread:num_elems, k) = AqTemp(i:num_elem_spread:num_elems, k);

            qTemp(i:num_elem_spread:num_elems, k) = 0.0;
        end
    end
end
