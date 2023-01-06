function [vars] = run_model(fn_input)
    %% Read namelist, and create directory to save output.
    addpath('util')    
    fprintf('Reading namelist %s...\n', fn_input)
    vars = read_namelist(fn_input);
    mkdir(sprintf('%s/%s', vars.base_dir, vars.exp_name));
    copyfile(fn_input, sprintf('%s/%s/namelist_%s.txt', vars.base_dir, vars.exp_name, vars.exp_name));
    %% Set up the parameters of the model. 
    % In RCE, kappa (as defined in Khairoutdinov and Emanuel, 2018), must
    % be equal to 1
    vars.kappa = 1; 
    vars.ep = 0.8;                              % precipitation efficiency
    
    %%% Time stepping settings
    vars.t = linspace(0, vars.T, vars.n_steps); % time stamps for each step
    vars.dt = vars.t(2) - vars.t(1);            % time-stepping resolution

    %%% Horizontal structure settings
    % Zonal resolution is only used for plotting.
    % Meridional resolution is the numerical resolution in the y-direction.
    vars.x = linspace(-pi, pi, vars.Nx);        % x dimension (for plots)
    vars.y = linspace(-3*pi, 3*pi, vars.Ny);    % y dimension
    vars.dy = vars.y(2) - vars.y(1);            % grid-spacing in y

    %%% Vertical structures settings
    % Troposphere vertical resolution is only used for plotting.
    vars.tr_p = flip(linspace(vars.p_t, vars.p_s, vars.Np))*100;    % troposphere pressure (Pa)
    vars.tr_p_nd = vars.tr_p ./ (vars.tr_p(1) - vars.tr_p(end));    % troposphere pressure (non-dimensional)
    
    % Generate a moist adiabat starting at the surface temperature.
    % This can be done by accumulating the difference from the mean temperature.
    e_sat = @(T) 6.1094 * exp(17.625 * (T - 273) ./ (T - 273 + 243.04)) * 100;  % saturation vapor pressure (Pa)
    q_sat = @(T, p) 0.622 * e_sat(T) ./ (p - e_sat(T));                         % saturation specific humidity
    gamma_w = @(T, p) (1 / p) * (vars.Rd * T + vars.Lv * q_sat(T, p)) / (vars.cp_d + (vars.Lv.^2 * q_sat(T, p) * 0.622) / (vars.Rd * T.^2));
    vars.T_m = zeros(numel(vars.tr_p), 1);                                      % temperature (K) in troposphere
    vars.T_m(1) = vars.T_bl;
    for i = 1:(numel(vars.tr_p)-1) 
        vars.T_m(i+1) = vars.T_m(i) + gamma_w(vars.T_m(i), vars.tr_p(i)) * (vars.tr_p(i+1) - vars.tr_p(i));
    end
    vars.T_mn = -trapz(vars.tr_p_nd, vars.T_m);                                 % vertically averaged temperature (K)
    vars.Tv_m = vars.T_m ./ (1 - e_sat(vars.T_m) ./ vars.tr_p.' * (1 - 0.622)); % virtual temperature
    vars.Tv_mn = trapz(vars.tr_p, vars.Tv_m) ./ (vars.tr_p(end)- vars.tr_p(1)); % vertically averaged virtual temperature (K)
    vars.nuT_p = (vars.T_m - vars.T_mn) / (vars.T_bl - vars.T_mn);
    [~, vars.nuT_yp] = meshgrid(vars.y, vars.nuT_p);
    % Using hypsometric equation, relate to the thickness (use an adjusted vertically averaged
    % temperature to get the tropopause height exactly equal to the specified H).
    adj_T_mn = vars.H * 1000 / log(vars.tr_p_nd(1) ./ vars.tr_p_nd(end)) * vars.g / vars.Rd;     % K
    vars.tr_z = vars.Rd * adj_T_mn / vars.g * log(vars.tr_p_nd(1) ./ vars.tr_p_nd) / 1000;       % km

    % Stratosphere vertical resolution is the numerical resolution in stratosphere.
    % Uses high resolution near the tropopause, and low resolution in the upper stratosphere.
    vars.Nz = vars.Nz_lower + vars.Nz_upper;                        % stratosphere vertical resolution (log-pressure)
    vars.st_Tmn = 240;    
    vars.st_z_nd_lower = linspace(1, 1 + vars.res_edge, vars.Nz_lower + 1);
    vars.st_z_nd_lower = vars.st_z_nd_lower(1:end-1);
    vars.st_p_lower = exp(-(vars.st_z_nd_lower-1) * vars.H/vars.H_s).*vars.p_t*100;
    vars.st_z_lower = vars.H + vars.Rd * vars.st_Tmn / vars.g * log(vars.st_p_lower(1) ./ vars.st_p_lower) / 1000;
    vars.st_z_nd_upper = linspace(1 + vars.res_edge, vars.z_top, vars.Nz_upper);
    vars.st_p_upper = exp(-(vars.st_z_nd_upper-1) * vars.H/vars.H_s).*vars.p_t*100;
    vars.st_z_upper = vars.H + vars.Rd * vars.st_Tmn / vars.g * log(vars.st_p_lower(1) ./ vars.st_p_upper) / 1000;
    vars.st_z_nd = [vars.st_z_nd_lower, vars.st_z_nd_upper];        % log-pressure stratosphere coordinates
    vars.st_p = exp(-(vars.st_z_nd-1) * vars.H/vars.H_s).*vars.p_t*100;
    vars.st_rho_nd = exp(vars.H / vars.H_s * (1 - vars.st_z_nd));   % non-dimensional density in stratosphere
    [~, vars.st_rho_nd_yz] = meshgrid(vars.y, vars.st_rho_nd);
    vars.st_z = vars.tr_z(end) + vars.Rd * vars.st_Tmn / vars.g * log(vars.st_p(1) ./ vars.st_p) / 1000; % km
    [vars.Y_YZ, vars.Z_YZ] = meshgrid(vars.y, vars.st_z_nd);        % stratosphere grid

    %%% Calculate dimensional variables from various non-dimensional variables in order
    %%% to convert back to dimensional units.
    vars.sd = vars.cp_d * log(vars.T_m / vars.T_m(1)).' - vars.Rd * log(vars.tr_p / vars.tr_p(1));  % dry entropy (J / K)
    vars.B = vars.H_s / vars.H * (vars.p_s - vars.p_t) ./ vars.p_t;
    vars.dsd_dz = (vars.sd(end) - vars.sd(1)) / (vars.H * 1000);                                    % dry entropy gradient (J / K m)
    vars.Ly = (1.7 * vars.dsd_dz * (vars.T_bl - vars.T_mn) * (vars.H * 1000) * (1 - vars.ep) / vars.beta^2)^0.25;   % deformation radius
    vars.B_w_to_omega = vars.T_m.' * vars.Rd / (vars.g * vars.H * 1000) * (vars.p_s - vars.p_t) ./ (vars.tr_p / 100);
    vars.Tstrat_mn = vars.H_s * 1000 * vars.g / vars.Rd;

    % Initialize sponge layer. They are active at the left, right, and top
    sponge_tr_P = vars.r_sponge*max(0, ((vars.y - vars.y_sponge) ./ (vars.y(end) - vars.y_sponge)).^4);
    sponge_tr_P(vars.y < vars.y_sponge) = 0;
    sponge_tr_N = vars.r_sponge*max(0, ((vars.y + vars.y_sponge) ./ (vars.y(1) + vars.y_sponge)).^4);
    sponge_tr_N(vars.y > -vars.y_sponge) = 0;
    vars.sponge_tr = sponge_tr_P + sponge_tr_N;
    
    sponge_tr_P = vars.r_sponge*max(0, ((vars.Y_YZ - vars.y_sponge) ./ (vars.Y_YZ(end) - vars.y_sponge)).^4);
    sponge_tr_P(vars.Y_YZ < vars.y_sponge) = 0;
    sponge_tr_N = vars.r_sponge*max(0, ((vars.Y_YZ + vars.y_sponge) ./ (vars.Y_YZ(1) + vars.y_sponge)).^4);
    sponge_tr_N(vars.Y_YZ > -vars.y_sponge) = 0;
    vars.sponge_st = sponge_tr_P + sponge_tr_N;
    vars.sponge_st((abs(vars.Z_YZ) < vars.z_sponge) & (abs(vars.Y_YZ) < 3)) = 0;
    vars.sponge_st = vars.sponge_st + vars.r_sponge*max(0, ((vars.Z_YZ - vars.z_sponge) ./ (vars.st_z_nd(end)-vars.z_sponge)).^3);

    %% Initialize vertical profile of mean wind in the stratosphere.
    [~, min_idx] = min(abs(vars.st_z - vars.H_uS_const));
    if strcmpi(vars.uS_shape, 'jump')
        vars.mean_Us_yz = vars.Us * ones(size(vars.Z_YZ));
    elseif strcmpi(vars.uS_shape, 'linear')
        if vars.Us < 0
            vars.mean_Us_yz = max(vars.Us, vars.Us *(vars.Z_YZ - 1) / ((vars.st_z_nd(min_idx) - 1)));
        else
            vars.mean_Us_yz = min(vars.Us, vars.Us * (vars.Z_YZ - 1) / ((vars.st_z_nd(min_idx) - 1)));
        end
        vars.mean_dTdy = zeros(size(vars.mean_Us_yz));
    elseif strcmpi(vars.uS_shape, 'qbo')
        vars.reduc_fac = 1 - exp(-(vars.Z_YZ - 1).^2 / vars.c_qbo);
        vars.mean_Us_yz = vars.Us * sin(vars.a_qbo * (vars.Z_YZ-1)) .* exp(-vars.b_qbo*vars.Y_YZ.^2) .* vars.reduc_fac;
        vars.mean_dUsdz = vars.a_qbo * vars.Us * cos(vars.a_qbo * (vars.Z_YZ-1)) .* exp(-vars.b_qbo*vars.Y_YZ.^2) .* vars.reduc_fac + ...
                          2 * vars.Us * (vars.Z_YZ-1).^3 .* exp(-(vars.Z_YZ - 1).^4 / vars.c_qbo) .* sin(vars.a_qbo .* (vars.Z_YZ - 1)) / vars.c_qbo .* exp(-vars.b_qbo*vars.Y_YZ.^2);
        vars.mean_T_yz = vars.mean_dUsdz * vars.H_s / vars.H / (2 * vars.b_qbo);
        vars.mean_dTdy = dfdy_2d(vars.mean_T_yz, vars, 2);
    else
        error('Not a valid stratosphere wind profile.')
    end

    %% Initialize the domain, either from an existing file or rigid-lid solution.
    if vars.init_from_file
        fprintf('Using initialization file %s...\n', vars.init_file);
        % TODO: Interpolate from the initialization file. May throw an
        % error right now if the current grid is not the same as the initizliation grid.
        init_mat = load(vars.init_file);
        state_matrix = init_mat.state_matrix;
        vars.sig_rigid = init_mat.sig_rigid;
        copyfile(vars.init_file, sprintf('%s/%s/n%dk%d_%s_init.mat', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name));
    else
        fprintf('Generating initialization file from rigid-lid solution...\n')
        vars.init_file = sprintf('%s/%s/n%dk%d_%s_init.mat', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name);
        vars.sig_rigid = equa_init(vars, 1);
        init_mat = load(vars.init_file);
        state_matrix = init_mat.state_matrix;
    end
    vars.output_file = sprintf('%s/%s/n%dk%d_%s_output.mat', vars.base_dir, vars.exp_name, vars.nm, vars.km, vars.exp_name);
    fprintf('Will be saving to output file %s...\n', vars.output_file);

    %%% Numerical computation.
    % Linear solutions will require rescaling.
    state_step = zeros(vars.n_steps / vars.save_state_step, size(state_matrix, 1), vars.Ny);
    rescale_step = zeros(vars.n_steps / vars.save_state_step, 1);
    vars.Cm = flip(dlmread('redblue_256_rgb.txt') / 255);   % colormap

    tic
    for it = 1:vars.n_steps
        % Print time elapsed every time a save step is reached.
        if mod(it, vars.save_state_step) == 0
            elapsedTime = toc;
            fprintf('Elapsed simulation time: %0.2f days; elapsed clock time: %2.2f minutes\n', ...
                    it / vars.n_steps * vars.T * vars.a / (vars.beta .* vars.Ly.^2) / 86400, ...
                    elapsedTime / 60);
        end
    
        % RK4 time-stepping.
        k1 = vars.dt * dState_dt(state_matrix, vars);      
        k2 = vars.dt * dState_dt(state_matrix + k1 / 2, vars);
        k3 = vars.dt * dState_dt(state_matrix + k2 / 2, vars);
        k4 = vars.dt * dState_dt(state_matrix + k3, vars);
        dstate_dt = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        state_matrix = state_matrix + dstate_dt;
        
        [uBT_y, uBC_y, vBT_y, vBC_y, ~, smT_y, uS_yz, vS_yz, phiS_yz, etaTp_y, ~] = matrix_to_state(state_matrix, vars);
        
        if vars.plot_steps && (mod(it, 100) == 0 || it == 1)
            % If plotting is turned on, plots the domain zonal velocity and
            % vertical velocity at the current time step.
            uBC_yp = baroclinic_vertical_mode(uBC_y, vars);
            uBT_yp = barotropic_vertical_mode(uBT_y, vars);
            
            omegaTp_y = 1i * vars.km * uBT_y + dfdy(vBT_y, vars, 2);
            wTp_y = -omegaTp_y * vars.B + 1i * vars.km * vars.mean_Us_yz(1, :) .* etaTp_y;
            [wTp_yz, ~] = meshgrid(wTp_y, vars.st_z_nd);
            wSt_yz = wTp_yz - cumtrapz(vars.st_z_nd, 1i * vars.km * uS_yz + dfdy_2d(vS_yz, vars, 2));

            u_ct = vertcat(uBT_yp + uBC_yp, uS_yz);
            z_ct = [vars.tr_z, vars.st_z];
            subplot(1, 2, 1);
            hold off; pcolor(vars.y, z_ct, real(u_ct)); shading interp; colorbar
            cmax = max(abs(u_ct), [], 'all'); caxis([-cmax, cmax]);  ylim([0, 50]);
            xlim([vars.y(1), vars.y(end)]); set(gca, 'YDir', 'Normal')
            xlabel('y (non-dimensional)'); ylabel('z (non-dimensional)');
            title(sprintf('Time: %0.2f, Time step %d', vars.t(it), it));
            
            subplot(1, 2, 2);
            hold off; pcolor(vars.y, vars.st_z, real(wSt_yz)); shading interp; colorbar;
            ylim([vars.st_z(1), 100]); cmax = max(abs(wSt_yz), [], 'all'); caxis([-cmax, cmax]); 
            xlim([vars.y(1), vars.y(end)]); set(gca, 'YDir', 'Normal')
            xlabel('y (non-dimensional)'); ylabel('z (non-dimensional)');
            colormap(vars.Cm); title('Stratospheric w');
            pause(0.01);
        end
        
        if mod(it, vars.save_state_step) == 0
            istep = it / vars.save_state_step;
            state_step(istep, :, :) = state_matrix;
            
            % Rescale everytime the mean domain energy is above the
            % rescale threshold.
            if mean(abs(phiS_yz).^2, 'all') > vars.rescale_threshold
                state_matrix = state_matrix / vars.rescale_threshold;
                vars.D_st = vars.D_st / vars.rescale_threshold;
                rescale_step(istep) = 1;
            end
        end
    end

    % Save the simulation output, along with steps where the domain rescaled.
    save(vars.output_file, 'state_step', 'rescale_step', 'vars', '-v7.3');
    fprintf('Saved %s!\n', vars.output_file);

    % Analyze output.
    fprintf('Analyzing output...\n');
    equa_analyze(vars);
end
