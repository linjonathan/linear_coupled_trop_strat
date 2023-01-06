function [dFdz] = dfdz_2d(f, dz, l_bc)    
    % y-derivative (dimension 2) of a 2-D matrix.
    % l_bc = 0: dfdy(1) = 0; symmetric    
    % l_bc = 1: f(1) = 0; antisymmetric
    
    dFdz = zeros(size(f));
    coeffs = [1/12, -2/3, 2/3, -1/12] / dz;
    dFdz(3:end-2, :) = coeffs(1) * f(1:end-4, :) + coeffs(2) * f(2:end-3, :) + ...
                       coeffs(3) * f(4:end-1, :) + coeffs(4) * f(5:end, :);    

    % Use boundary conditions for edges.
    if l_bc == 0
        % Symmetric; f(1) = f(-1) and f(2) = f(-2)
        dFdz(1, :) = 0;
        dFdz(2, :) = coeffs(1) * f(2, :) + coeffs(2) * f(1, :) + ...
                  coeffs(3) * f(3, :) + coeffs(4) * f(4, :);
    elseif l_bc == 1
        % Antisymmetric; f(1) = -f(-1) and f(2) = -f(-2)
        dFdz(1, :) = 2*coeffs(3) * f(2, :) + 2*coeffs(4) * f(3, :);
        dFdz(2, :) = -coeffs(1) * f(2, :) + coeffs(2) * f(1, :) + ...
                  coeffs(3) * f(3, :) + coeffs(4) * f(4, :);
    else
        %dFdz(1, :) = (1/3*f(4,:) - 3/2*f(3,:) + 3*f(2, :) - 11/6*f(1, :)) ./ dz;
        dFdz(1, :) = (-1/2*f(3,:) + 2*f(2, :) - 3/2*f(1, :)) ./ dz;
        dFdz(2, :) = (f(3, :) - f(1, :)) ./ (2 * dz);
    end
    dFdz(end-1, :) = (f(end, :) - f(end-2, :)) / (2 * dz);    
    dFdz(end, :) = (f(end, :) - f(end-1, :)) ./ dz;    
end