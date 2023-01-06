function [] = equa_analyze(vars)
    % Generate various figures and analysis on the linear solution.
    out_mat = load(vars.output_file);
    state_step = out_mat.state_step;
    rescale_step = out_mat.rescale_step;
    vars = out_mat.vars;
    
    its = 1:vars.n_steps;
    t_step = vars.t(mod(its, vars.save_state_step) == 0);
    rescale_steps = find(rescale_step==1);
    sigs = zeros(1, numel(vars.y)); sigs(:) = NaN;
    no_sponge_idxs = find(vars.sponge_tr(1, :)  < (vars.r_sponge / 10));
    T_nd_start = vars.T_analyze_start;
    T_nd_end = vars.T_analyze_end;

    for y_idx = no_sponge_idxs
        % Use s to estimate growth rate.
        s_log_r = log(abs(state_step(:, 5, y_idx)).');
        s_thta_all = atan2(imag(state_step(:, 5, y_idx)), real(state_step(:, 5, y_idx)));
        
        t_start = find(t_step >= T_nd_start, 1);
        t_end = find(t_step >= T_nd_end, 1);
        mask = (rescale_steps >= t_start) & (rescale_steps <= t_end);
        t_rescale_valid = rescale_steps(mask);
        for ts_idx = 1:numel(t_rescale_valid)
            ts = t_rescale_valid(ts_idx);
            if ts < numel(s_log_r)
                s_log_r((ts+1):end) = s_log_r((ts+1):end) - s_log_r(ts+1) ...
                    + s_log_r(ts) + (s_log_r(ts) - s_log_r(ts-1));
            end
        end
        s_log_r = s_log_r((t_step >= t_step(t_start)) & (t_step <= t_step(t_end)));
        t_mask = (t_step >= t_step(t_start)) & (t_step <= t_step(t_end));
        t_var = t_step(t_mask);
        linear_reg = polyfit(t_var, s_log_r, 1);
        sig_r = linear_reg(1);


        s_thta = s_thta_all(t_mask);
        s_diff = s_thta(2:end) - s_thta(1:end-1);
        pass_idxs = find(abs(s_diff) > (pi/2));
        s_thta_corr = s_thta;
        first_idx = find(t_mask, 1);
        for ip = 1:numel(pass_idxs)
            p_idx = pass_idxs(ip);
            if p_idx == 1
                s_thta_corr((p_idx+1):end) = s_thta_corr((p_idx+1):end) - s_thta_corr(p_idx+1) ...
                         + s_thta_corr(p_idx) + (s_thta_corr(p_idx) - s_thta_all(first_idx-1));

            else
                s_thta_corr((p_idx+1):end) = s_thta_corr((p_idx+1):end) - s_thta_corr(p_idx+1) ...
                         + s_thta_corr(p_idx) + (s_thta_corr(p_idx) - s_thta_corr(p_idx-1));
            end
        end
        linear_reg = polyfit(t_var.', s_thta_corr, 1);
        sig_i = abs(linear_reg(1));
        sigs(y_idx) = sig_r - 1i * sig_i;
    end

    figure;
    subplot(2, 1, 1);  plot(vars.y, real(sigs)); hold on;
    plot(vars.y, ones(numel(vars.y), 1) * real(vars.sig_rigid), 'k--');
    xlabel('y (non-dimensional)'); ylabel('Magnitude')
    set(gca, 'Fontsize', 13); grid on
    xlim([-5, 5]);
    legend(sprintf('S = %d', vars.S), 'Rigid Lid');title('Real(\sigma)');

    subplot(2, 1, 2); plot(vars.y, imag(sigs)); hold on;
    plot(vars.y, ones(numel(vars.y), 1) * imag(vars.sig_rigid), 'k--');
    xlabel('y (non-dimensional)'); ylabel('Magnitude')
    set(gca, 'Fontsize', 13); grid on
    xlim([-5, 5]);
    legend(sprintf('S = %d', vars.S), 'Rigid Lid'); title('Imag(\sigma)'); 
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_sig.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');

    % Check solutions. Validation of the modes.        
    [~, e_idx] = min(abs(vars.y));
    sig = sigs(e_idx);
    fac = exp(imag(sig)*1i*t_step);
    tc_start = find(t_step >= T_nd_start, 1); tc_end = find(t_step >= T_nd_end, 1);
    state_end = squeeze(state_step(tc_start:tc_end, :, :));
    state_centered = zeros(size(state_end));
    for it = 1:(tc_end - tc_start + 1)
        state_center = state_end(it, :, :) / fac(it);
        state_centered(it, :, :) = state_center;  
    end
    for it = 1:(tc_end - tc_start + 1)
        state_centered(it, :, :) = state_centered(it, :, :) ./ abs(squeeze(state_centered(it, 2, e_idx)));
    end
    y_structs = squeeze(mean(state_centered, 1));        
    y_structs(isnan(y_structs)) = 0;
    [uBT_y, uBC_y, vBT_y, vBC_y, ~, ~, uS_yz, vS_yz, phiS_yz, etaTp_y, ~] = matrix_to_state(y_structs, vars);
    uBC_yp = baroclinic_vertical_mode(uBC_y, vars);
    vBC_yp = baroclinic_vertical_mode(vBC_y, vars);
    uBT_yp = barotropic_vertical_mode(uBT_y, vars);
    vBT_yp = barotropic_vertical_mode(vBT_y, vars);    

    omegaTp_y = 1i * vars.km * uBT_y + dfdy(vBT_y, vars, 2);
    wTp_y = -omegaTp_y * vars.B + 1i * vars.km * vars.mean_Us_yz(1, :) .* etaTp_y;
    [wTp_yz, ~] = meshgrid(wTp_y, vars.st_z_nd);    
    wSt_yz = wTp_yz - cumtrapz(vars.st_z_nd, 1i * vars.km * uS_yz + dfdy_2d(vS_yz, vars, 2));

    dstatedt = dState_dt(y_structs, vars);
    n_vars = 9;
    var_names = {'u_{bt}', 'u_{bc}', 'v_{bt}', 'v_{bc}', 's_{trop}', 'sm_{trop}', 'u_{strat}', 'v_{strat}', '\phi_{tp}'};

    figure('units','normalized','outerposition',[0.1 0.1 0.6 0.8]);
    for i_var = 1:n_vars
        if sign(imag(sig)) == -1
            str_sgn = '-';
        else
            str_sgn = '+';
        end    
        subplot(5, 2, i_var); hold on;
        if i_var <= 7
            plot(vars.y, real(dstatedt(i_var, :)), 'b', 'Linewidth', 1.5); 
            plot(vars.y, imag(dstatedt(i_var, :)), 'r', 'Linewidth', 1.5); 
            plot(vars.y, real(sig * y_structs(i_var, :)), 'g--', 'Linewidth', 0.5);
            plot(vars.y, imag(sig * y_structs(i_var, :)), 'g-.', 'Linewidth', 0.5);
        else
            idx_offset = 7; idx_var = vars.Nz;
            plot(vars.y, real(dstatedt(idx_offset + idx_var*(i_var-7), :)), 'b', 'Linewidth', 1.5); 
            plot(vars.y, imag(dstatedt(idx_offset + idx_var*(i_var-7), :)), 'r', 'Linewidth', 1.5); 
            plot(vars.y, real(sig * y_structs(idx_offset + idx_var*(i_var-7), :)), 'g--', 'Linewidth', 0.5);
            plot(vars.y, imag(sig * y_structs(idx_offset + idx_var*(i_var-7), :)), 'g-.', 'Linewidth', 0.5);
        end
        title(sprintf('%s: %f %s %fi', var_names{i_var}, real(sig), str_sgn, abs(imag(sig))));
    end
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_eigenmode_y.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');

    figure;
    pcolor(vars.y, vars.st_z, real(wSt_yz)); shading interp
    colorbar; cmax = max(abs(wSt_yz), [], 'all');
    caxis([-cmax, cmax]); colormap(vars.Cm);
    ylabel('Height (km)'); ylim([vars.H, 80]);
    xlabel('y (non-dimensional)');
    set(gca, 'Fontsize', 14);
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_w_yz.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png'); 

    figure('units','normalized','outerposition',[0.3 0.1 0.4 0.8]);
    for i_var = 1:n_vars
        subplot(5, 2, i_var);
        x = vars.x;
        y = vars.y;
        [X, ~] = meshgrid(x, y);
        if i_var <= 7
            var_mode = exp(1i*X) .* repmat(squeeze(y_structs(i_var, :)).', 1, vars.Nx);
        else
            idx_offset = 7; idx_var = vars.Nz;
            var_mode = exp(1i*X) .* repmat(squeeze(y_structs(idx_offset + idx_var*(i_var-7), :)).', 1, vars.Nx);
        end                     
        pcolor(x, y, real(var_mode)); shading interp; colorbar;
        cmax = max(abs(var_mode), [], 'all');    
        title(sprintf('%s', var_names{i_var}));
        ylim([-5, 5]); colormap(vars.Cm); caxis([-cmax, cmax]);
    end
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_eigenmode_xy.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');

    [X, ~] = meshgrid(vars.x, vars.y);
    s = exp(1i*X) .* repmat(squeeze(y_structs(5, :)).', 1, vars.Nx);
    u_bt = exp(1i*X) .* repmat(squeeze(y_structs(1, :)).', 1, vars.Nx);
    u_bc = exp(1i*X) .* repmat(squeeze(y_structs(2, :)).', 1, vars.Nx);
    v_bt = exp(1i*X) .* repmat(squeeze(y_structs(3, :)).', 1, vars.Nx);
    v_bc = exp(1i*X) .* repmat(squeeze(y_structs(4, :)).', 1, vars.Nx);
    w = -1i * vars.km * (u_bt + u_bc) - dfdy_2d(v_bt.' + v_bc.', vars, 2).';

    figure('units','normalized','outerposition',[0.0 0.0 0.4 0.5]);
    set(gcf,'color','w'); hold on;
    imagesc(x, vars.y, real(w));
    x_skip = 8; y_skip = 2;
    smax = max(real(s), [], 'all'); smax_incr = smax / 10;
    contour(x, vars.y, real(s), smax_incr:smax_incr:(smax+smax_incr), 'LineColor', [50 50 50] ./ 255);
    contour(x, vars.y, real(s), -(smax+smax_incr):smax_incr:-smax_incr, 'LineColor', [50 50 50] ./ 255, 'LineStyle', '--');
    U = real(u_bt + u_bc); V = real(v_bt + v_bc);
    quiver(x(1:x_skip:end), vars.y(1:y_skip:end), U(1:y_skip:end, 1:x_skip:end), V(1:y_skip:end, 1:x_skip:end), 'k'); % 'AutoScaleFactor', 0.3, 'MaxHeadSize', 0.5);
    xlim([x(1), x(end)]); ylim([-5, 5]);
    xlabel('x (non-dimensional)'); ylabel('y (non-dimensional)');
    cmax = max(real(w), [], 'all');
    colormap(vars.Cm); caxis([-cmax, cmax]);
    set(gca, 'Fontsize', 15);
    title('s (contour), u (arrows), w (color)');
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_troposphere_summary_xy.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');    

    % Stratosphere eigenmode
    z_lvls = [vars.H, vars.H+1, vars.H+2, 20, 25, 30];
    [X, ~] = meshgrid(vars.x, vars.y);
    for z_lvl = z_lvls
        [~, s_idx] = min(abs(vars.st_z - z_lvl));
        uS_xy = exp(1i*X) .* repmat(squeeze(y_structs(7 + s_idx - 1, :)).', 1, vars.Nx);
        vS_xy = exp(1i*X) .* repmat(squeeze(y_structs(idx_offset + idx_var + s_idx - 1, :)).', 1, vars.Nx);
        wSt_xy = exp(1i*X) .* repmat(squeeze(wSt_yz(s_idx, :)).', 1, vars.Nx); 
        mean_Us_xy = repmat(vars.mean_Us_yz(1, :).', 1, vars.Nx);
        tempS_xy = (-wSt_xy .* vars.S) ./ (sig + mean_Us_xy .* 1i * vars.km);
        
        figure('units','normalized','outerposition',[0.0 0.0 0.4 0.5]);
        set(gcf,'color','w'); hold on;
        imagesc(x, vars.y, real(wSt_xy));
        x_skip = 8; y_skip = 2;
        phimax = max(real(tempS_xy), [], 'all'); phi_incr = phimax / 10;
        contour(x, vars.y, real(tempS_xy), phi_incr:phi_incr:(phimax+phi_incr), 'LineColor', [50 50 50] ./ 255);
        contour(x, vars.y, real(tempS_xy), -(phimax+phi_incr):phi_incr:-phi_incr, 'LineColor', [50 50 50] ./ 255, 'LineStyle', '--');
        U = real(uS_xy); V = real(vS_xy);
        quiver(x(1:x_skip:end), vars.y(1:y_skip:end), U(1:y_skip:end, 1:x_skip:end), V(1:y_skip:end, 1:x_skip:end), 'k'); % 'AutoScaleFactor', 0.3, 'MaxHeadSize', 0.5);
        xlim([x(1), x(end)]); ylim([-5, 5]);
        xlabel('x (non-dimensional)'); ylabel('y (non-dimensional)');
        cmax = max(abs(wSt_xy), [], 'all');
        colormap(vars.Cm); caxis([-cmax, cmax]);
        set(gca, 'Fontsize', 15);
        title('T (contour), u (arrows), w (color)');
        saveas(gcf, sprintf('%s/%s/n%dk%d_%s_stratosphere_%dkm_summary_xy.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name, z_lvl), 'png');
    end
    
    figure('units','normalized','outerposition',[0.0 0.0 0.6 0.8]); set(gcf,'color','w');
    [X_tr, ~] = meshgrid(vars.x, vars.tr_p_nd);
    [X_st, ~] = meshgrid(vars.x, vars.st_z_nd);
    u_bt = exp(1i*X_tr) .* repmat(squeeze(uBT_yp(:, e_idx)), 1, vars.Nx);
    u_bc = exp(1i*X_tr) .* repmat(squeeze(uBC_yp(:, e_idx)), 1, vars.Nx);
    u_s = exp(1i*X_st) .* repmat(squeeze(uS_yz(:, e_idx)), 1, vars.Nx);
    v_bt = exp(1i*X_tr) .* repmat(squeeze(vBT_yp(:, e_idx)), 1, vars.Nx);
    v_bc = exp(1i*X_tr) .* repmat(squeeze(vBC_yp(:, e_idx)), 1, vars.Nx);
    v_s = exp(1i*X_st) .* repmat(squeeze(vS_yz(:, e_idx)), 1, vars.Nx);
        
    subplot(2, 2, 1);
    pcolor(vars.x, horzcat(vars.tr_z, vars.st_z), real(vertcat(u_bt + u_bc, u_s))); shading interp
    colorbar; cmax = max(real(u_bt + u_bc), [], 'all');
    caxis([-cmax, cmax]); colormap(vars.Cm);
    ylabel('Height (km'); ylim([0, 50]);
    xlabel('x (non-dimensional)'); 
    set(gca, 'Fontsize', 14); title('Zonal Velocity on Equator'); 
    subplot(2, 2, 2);
    pcolor(vars.x, horzcat(vars.tr_z, vars.st_z), real(vertcat(v_bt + v_bc, v_s))); shading interp
    colorbar; cmax = max(real(v_bt + v_bc), [], 'all');
    caxis([-cmax, cmax]); colormap(vars.Cm);
    ylabel('Height (km'); ylim([0, 50]);
    xlabel('x (non-dimensional)'); 
    set(gca, 'Fontsize', 14); title('Meridional Velocity on Equator'); 
        
    [~, e_idx] = min(abs(vars.y - 1.5));
    u_bt = exp(1i*X_tr) .* repmat(squeeze(uBT_yp(:, e_idx)), 1, vars.Nx);
    u_bc = exp(1i*X_tr) .* repmat(squeeze(uBC_yp(:, e_idx)), 1, vars.Nx);
    u_s = exp(1i*X_st) .* repmat(squeeze(uS_yz(:, e_idx)), 1, vars.Nx);
    v_bt = exp(1i*X_tr) .* repmat(squeeze(vBT_yp(:, e_idx)), 1, vars.Nx);
    v_bc = exp(1i*X_tr) .* repmat(squeeze(vBC_yp(:, e_idx)), 1, vars.Nx);
    v_s = exp(1i*X_st) .* repmat(squeeze(vS_yz(:, e_idx)), 1, vars.Nx);        
    
    subplot(2, 2, 3);
    pcolor(x, horzcat(vars.tr_z, vars.st_z), real(vertcat(u_bt + u_bc, u_s))); shading interp
    colorbar; cmax = max(real(u_bt + u_bc), [], 'all');
    caxis([-cmax, cmax]); colormap(vars.Cm);
    ylabel('Height (km'); ylim([0, 50]);
    xlabel('x (non-dimensional)'); 
    set(gca, 'Fontsize', 14); title('Zonal Velocity on y = 1.5'); 
    subplot(2, 2, 4);
    pcolor(x, horzcat(vars.tr_z, vars.st_z), real(vertcat(v_bt + v_bc, v_s))); shading interp
    colorbar; cmax = max(real(v_bt + v_bc), [], 'all');
    caxis([-cmax, cmax]); colormap(vars.Cm);
    ylabel('Height (km'); ylim([0, 50]);
    xlabel('x (non-dimensional)'); 
    set(gca, 'Fontsize', 14); title('Meridional Velocity on y = 1.5');
    saveas(gcf, sprintf('%s/%s/n%dk%d_%s_wind_xz.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');    

    figure;
    if abs(vars.Us) > 1e-2
        plot(vars.mean_Us_yz(:, 64), vars.st_z);
        ylabel('Height (km)');
        xlabel('Mean Wind (n.d.)');
        xlim([-abs(vars.Us)*1.1, abs(vars.Us)*1.1])
        set(gca, 'Fontsize', 13);  grid on;
        saveas(gcf, sprintf('%s/%s/n%dk%d_%s_strat_wind.png', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), 'png');
    end
  
    % Calculate wave energy flux
    [X, ~] = meshgrid(vars.x, vars.y);
    eFlux_yz = zeros(size(uS_yz));
    for iz = 1:size(uS_yz, 1)
        [~, wS_xy] = meshgrid(vars.x, wSt_yz(iz, :));
        [~, phiS_xy] = meshgrid(vars.x, phiS_yz(iz, :));
        wS_xy = wS_xy .* exp(1i * vars.km .* X) / norm(uBC_y, 2);
        phiS_xy = phiS_xy .* exp(1i * vars.km .* X) / norm(uBC_y, 2);
        eFlux_yz(iz, :) = mean(real(phiS_xy) .* real(wS_xy), 2);
    end
    eFlux_total = trapz(vars.y, trapz(vars.st_z_nd, eFlux_yz));

    save(sprintf('%s/%s/n%dk%d_%s_eFlux.mat', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name), ...
         'vars', 'eFlux_yz', 'eFlux_total', 'y_structs', 'sig');
end