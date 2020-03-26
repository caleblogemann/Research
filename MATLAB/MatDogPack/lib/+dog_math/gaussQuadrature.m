function [integral] = gaussQuadrature(integrandFunction, quad_order, lower_bound, upper_bound)
    [quad_pts, quad_wgts] = dog_math.getGaussQuadPtsAndWgts(quad_order);
    interval_width = (upper_bound - lower_bound);
    interval_center = 0.5*(upper_bound + lower_bound);
    integral = 0;
    for k = 1:quad_order
        % x = (b - a)/2*xi + (b + a)/2
        integrand_point = 0.5*interval_width*quad_pts(k) + interval_center;
        integral = integral + quad_wgts(k)*integrandFunction(integrand_point);
    end
end
