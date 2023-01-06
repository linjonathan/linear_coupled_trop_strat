function [dFdz] = dfdz_2d_stagger(f_yz, vars)
    dz_lower = vars.st_z_nd_lower(2) - vars.st_z_nd_lower(1);
    dz_upper = vars.st_z_nd_upper(2) - vars.st_z_nd_upper(1);
    dFdz_upper = dfdz_2d(f_yz((vars.Nz_lower + 1):end, :), dz_upper, 2);
    dFdz_lower = dfdz_2d(f_yz(1:vars.Nz_lower, :), dz_lower, 2);
    dFdz = vertcat(dFdz_lower, dFdz_upper);
end