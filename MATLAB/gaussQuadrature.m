function [integral] = gaussQuadrature(quad_order, integrandFunction)
    % assume integral between -1 and 1

    [quad_pts, quad_wgts] = getGaussQuadPtsAndWgts(quad_order);
    integral = 0;
    for k = 1:quad_order;
        integral = integral + quad_wgts(k)*integrandFunction(quad_pts(k));
    end
end
