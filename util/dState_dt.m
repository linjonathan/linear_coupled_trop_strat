function [dstate_dt] = dState_dt(state_matrix, vars)
    % Numerically integrates the coupled troposphere-stratosphere equations.
    % TROPOSPHERE EQUATIONS:
    %   du/dt = ds/dx + yv
    %   dv/dt = delta(ds/dy - yu)
    %   du/dx + dv/dy + w = 0
    %   ds/dt = C s_m - w + s_m - chi s - alpha u
    %   gamma * ds_m/dt = - DS - alpha u + kappa C s_m - G w + d d^2 s_m / dx^2
    % STRATOSPHERE EQUATIONS:
    %   du/dt = -dphi/dx + yv
    %   dv/dt = delta(-dphidy - yu)
    %   du/dx + dv/dy + 1 / S * d^2/dz^2 d/dt phi = 0;
    [uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y, uS_yz, vS_yz, phiS_yz, etaTp_y, qTp_y] = matrix_to_state(state_matrix, vars);
    
    % Tropopause variables.
    phiTp_y = phiS_yz(1, :);                                       % tropopause geopotential
    bcTp_y = vars.nuT_yp(end, :);                                  % baroclinic mode tropopause value
    omegaTp_y = 1i * vars.km * uBT_y + dfdy(vBT_y, vars, 2);

    % Use the vertical velocity matching condition across the tropopause.
    % If there is no jump in wind a the tropopause, enforce zero displacement.
    if sum(abs(vars.mean_Us_yz(1, :))) < 1e-5
        etaTp_y(:) = 0;
    end
    wTp_y = -omegaTp_y * vars.B + 1i * vars.km * vars.mean_Us_yz(1, :) .* etaTp_y;

    % Integrate mass continuity from the tropopause upwards.
    % Note, w should be near zero at the upper boundary.
    [wTp_yz, ~] = meshgrid(wTp_y, vars.st_z_nd);
    wSt_yz = wTp_yz - cumtrapz(vars.st_z_nd, 1i * vars.km * uS_yz + dfdy_2d(vS_yz, vars, 2));

    % Geopotential in the stratosphere
    dphiSdz = dfdz_2d_stagger(phiS_yz, vars);
    if vars.z_cirrus >= vars.H
        idx_tp = find(vars.st_z >= vars.z_cirrus, 1);
        wCirrus_y = wSt_yz(idx_tp, :);
    else
        idx_tp = find(vars.tr_z >= vars.z_cirrus, 1);
        uBC_yp = baroclinic_vertical_mode(uBC_y, vars);
        vBC_yp = baroclinic_vertical_mode(vBC_y, vars);
        uBT_yp = barotropic_vertical_mode(uBT_y, vars);
        vBT_yp = barotropic_vertical_mode(vBT_y, vars);
        uT_yp = uBC_yp + uBT_yp;
        vT_yp = vBC_yp + vBT_yp;
        [~, BWO] = meshgrid(vars.y, vars.B_w_to_omega);
        [~, BWW] = meshgrid(vars.y, vars.Tstrat_mn ./ vars.T_m);
        wT_yp = -BWW .* BWO .* cumtrapz(vars.tr_p_nd, -1i * vars.km * uT_yp - dfdy(vT_yp.', vars, 2).');
        wCirrus_y = wT_yp(idx_tp, :);
    end
    u_advect_qTp = vars.u_advect;

    % Boundary layer variables.
    % Sum of barotropic and baroclinic since baroclinic mode is equal to one in boundary layer.
    uBL_y = uBT_y + uBC_y;
    vBL_y = vBT_y + vBC_y;
    
    % Linearizing about mean easterlies leads to factor 2 in zonal
    % friction.
	F_fac = 2;

    % Troposphere equations
    phiBL_y = phiTp_y - (1 - bcTp_y) .* sT_y;
    phiBT_y = phiBL_y + sT_y;
    duBTdt = -1i * vars.km * phiBT_y + vars.y .* vBT_y - vars.sponge_tr .* uBT_y ...
	     - F_fac * vars.F * uBL_y + vars.eta * dfdy(dfdy(uBT_y, vars, 2), vars, 2);
    duBCdt = 1i * vars.km * sT_y + vars.y .* vBC_y - vars.sponge_tr .* uBC_y ...
	     - F_fac * vars.F * uBL_y + vars.eta * dfdy(dfdy(uBC_y, vars, 2), vars, 2);
    dvBTdt = vars.delta * (-dfdy(phiBT_y, vars, 2) - vars.y .* uBT_y) - vars.F * vBL_y ...
             - vars.sponge_tr .* vBT_y + vars.eta * dfdy(dfdy(vBT_y, vars, 2), vars, 2);
    dvBCdt = vars.delta * (dfdy(sT_y, vars, 2) - vars.y .* uBC_y) - vars.F * vBL_y ...
             - vars.sponge_tr .* vBC_y + vars.eta * dfdy(dfdy(vBC_y, vars, 2), vars, 2);
    wT_y = -(1i * vars.km * uBL_y + dfdy_2d(vBL_y, vars, 2));
    dsTdt = (1+vars.C) * smT_y + vars.Cr * qTp_y - wT_y - vars.chi * sT_y - vars.alpha * uBL_y - vars.sponge_tr .* sT_y + ...
             vars.eta * dfdy(dfdy(sT_y, vars, 2), vars, 2);
    dsmTdt = 1 / vars.gamma * (-vars.D .* sT_y - vars.alpha .* uBL_y + vars.kappa * vars.C * smT_y + vars.Cr * qTp_y ...
             -vars.G * wT_y) - vars.sponge_tr .* smT_y + vars.eta * dfdy(dfdy(smT_y, vars, 2), vars, 2);
    dqTpdt = -1i * vars.km * u_advect_qTp .* qTp_y + vars.upsilon * wCirrus_y;
    
     % Stratosphere equations
    duSdt = -1i * vars.km * phiS_yz ...
            + vars.Y_YZ .* vS_yz ...
            - vars.mean_Us_yz .* 1i * vars.km .* uS_yz ...
            - vars.Gamma * wSt_yz .* dfdz_2d_stagger(vars.mean_Us_yz, vars) ...
            - vars.sponge_st .* uS_yz ...
	        + vars.eta * dfdy_2d(dfdy_2d(uS_yz, vars, 2), vars, 2);
    dvSdt = vars.delta * (-dfdy_2d(phiS_yz, vars, 2) - vars.Y_YZ .* uS_yz) ...
            - vars.mean_Us_yz .* 1i * vars.km .* vS_yz ...
	        - vars.sponge_st .* vS_yz  ...
	        + vars.eta * dfdy_2d(dfdy_2d(vS_yz, vars, 2), vars, 2);

    dphiSdt = - 1i * vars.km * flip(cumtrapz(flip(vars.st_z_nd), flip(vars.mean_Us_yz .* dphiSdz))) ...
              - flip(cumtrapz(flip(vars.st_z_nd), flip(wSt_yz) * vars.S)) ...
              - vars.xi * flip(cumtrapz(flip(vars.st_z_nd), flip(vars.mean_dTdy .* vS_yz))) ...
              - vars.q_rad * phiS_yz ...
	          - vars.sponge_st .* phiS_yz ...
	          + vars.eta * dfdy_2d(dfdy_2d(phiS_yz, vars, 2), vars, 2);

    % Tropopause equations
    detaTpdt = -omegaTp_y * vars.B;
    
    % Spectral damping to get rid of unwanted noise.
    l = [(-vars.Ny/2):((vars.Ny/2)-1)];
    spectral_damp_y = -tanh(abs(l / 2.5) - 10) / 2  + 0.5;
    [spec_damp_YZ, ~] = meshgrid(spectral_damp_y, vars.st_z_nd);
    IuBT_y = fftshift(fft(duBTdt)) .* spectral_damp_y;
    IuBC_y = fftshift(fft(duBCdt)) .* spectral_damp_y;
    IvBT_y = fftshift(fft(dvBTdt)) .* spectral_damp_y;
    IvBC_y = fftshift(fft(dvBCdt)) .* spectral_damp_y;
    IsT_y = fftshift(fft(dsTdt)) .* spectral_damp_y;
    IsmT_y = fftshift(fft(dsmTdt)) .* spectral_damp_y;    
    IuS_yz = fftshift(fft(duSdt, [], 2), 2) .* spec_damp_YZ;
    IvS_yz = fftshift(fft(dvSdt, [], 2), 2) .* spec_damp_YZ;
    IphiS_yz = fftshift(fft(dphiSdt, [], 2), 2) .* spec_damp_YZ;
    IetaTp_y = fftshift(fft(detaTpdt)) .* spectral_damp_y;
    IqTp_y = fftshift(fft(dqTpdt)) .* spectral_damp_y;
    duBTdt = ifft(ifftshift(IuBT_y));
    duBCdt = ifft(ifftshift(IuBC_y));
    dvBTdt = ifft(ifftshift(IvBT_y));
    dvBCdt = ifft(ifftshift(IvBC_y));
    dsTdt = ifft(ifftshift(IsT_y));
    dsmTdt = ifft(ifftshift(IsmT_y));
    duSdt = ifft(ifftshift(IuS_yz, 2), [], 2);
    dvSdt = ifft(ifftshift(IvS_yz, 2), [], 2);
    dphiSdt = ifft(ifftshift(IphiS_yz, 2), [], 2);
    detaTpdt = ifft(ifftshift(IetaTp_y));
    dqTpdt = ifft(ifftshift(IqTp_y));

    if sum(abs(vars.mean_Us_yz(1, :))) < 1e-5
        detaTpdt(:) = 0;
    end

    dstate_dt = state_to_matrix(duBTdt, duBCdt, dvBTdt, dvBCdt, dsTdt, ...
                                dsmTdt, duSdt, dvSdt, dphiSdt, detaTpdt, dqTpdt);
end
