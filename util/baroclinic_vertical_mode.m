function [fBC_yp] = baroclinic_vertical_mode(fBC_y, vars)
    [fBC_yp, ~] = meshgrid(fBC_y, vars.tr_p_nd); 
    fBC_yp = fBC_yp .* vars.nuT_yp;
end