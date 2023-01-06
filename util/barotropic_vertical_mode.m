function [fBT_yp] = barotropic_vertical_mode(fBT_y, vars)
    [fBT_yp, ~] = meshgrid(fBT_y, vars.tr_p_nd); 
end