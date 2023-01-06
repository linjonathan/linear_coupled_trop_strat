function [dFdy] = dfdy_2d(f, vars, l_bc)    
    % y-derivative (dimension 2) of a 2-D matrix.
    % l_bc = 0: dfdy(1) = 0; symmetric    
    % l_bc = 1: f(1) = 0; antisymmetric
    
    dFdy = zeros(size(f));
    coeffs = [1/12, -2/3, 2/3, -1/12] / vars.dy;
    dFdy(:, 3:end-2) = coeffs(1) * f(:, 1:end-4) + coeffs(2) * f(:, 2:end-3) + ...
                       coeffs(3) * f(:, 4:end-1) + coeffs(4) * f(:, 5:end);    

    % Use boundary conditions for edges.
    if l_bc == 0
        % Symmetric; f(1) = f(-1) and f(2) = f(-2)
        dFdy(:, 1) = 0;
        dFdy(:, 2) = coeffs(1) * f(:, 2) + coeffs(2) * f(:, 1) + ...
                  coeffs(3) * f(:, 3) + coeffs(4) * f(:, 4);
    elseif l_bc == 1
        % Antisymmetric; f(1) = -f(-1) and f(2) = -f(-2)
        dFdy(:, 1) = 2*coeffs(3) * f(:, 2) + 2*coeffs(4) * f(:, 3);
        dFdy(:, 2) = -coeffs(1) * f(:, 2) + coeffs(2) * f(:, 1) + ...
                  coeffs(3) * f(:, 3) + coeffs(4) * f(:, 4);
    else
        dFdy(:, 1) = (f(:, 2) - f(:, 1)) ./ vars.dy;
        dFdy(:, 2) = (f(:, 3) - f(:, 1)) ./ (2 * vars.dy);
    end
    dFdy(:, end-1) = (f(:, end) - f(:, end-2)) / (2 * vars.dy);    
    dFdy(:, end) = (f(:, end) - f(:, end-1)) ./ vars.dy;    
end