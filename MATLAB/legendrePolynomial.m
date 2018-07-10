function [phi] = legendrePolynomial(k, xi)
    if(k == 1)
        phi = 1;
    elseif (k==2)
        phi = sqrt(3)*xi;
    elseif (k==3)
        phi = sqrt(5)/2*(3*xi.^2 - 1);
    elseif (k==4)
        phi = sqrt(7)/2*(5*xi.^3 - 3*xi);
    else
        error('This order polynomial is not implemented');
    end
end
