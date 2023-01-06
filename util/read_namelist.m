function [vars] = read_namelist(fn_input)
    fileID = fopen(fn_input);
    tline = fgetl(fileID);          % header line
    tline = fgetl(fileID);          % header line
    vars.nm = sscanf(fgetl(fileID), 'nm=%f;');
    vars.km = sscanf(fgetl(fileID), 'km=%f;*');
    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.alpha = sscanf(fgetl(fileID), 'alpha=%f;');
    vars.chi = sscanf(fgetl(fileID), 'chi=%f;*');
    vars.C = sscanf(fgetl(fileID), 'C=%f;*');
    vars.Cr = sscanf(fgetl(fileID), 'Cr=%f;*');
    vars.gamma = sscanf(fgetl(fileID), 'gamma=%f;*');
    vars.D = sscanf(fgetl(fileID), 'D=%f;*');
    vars.d = sscanf(fgetl(fileID), 'd=%f;*');
    vars.G = sscanf(fgetl(fileID), 'G=%f;*'); 
    vars.delta = sscanf(fgetl(fileID), 'delta=%f;*');
    vars.F = sscanf(fgetl(fileID), 'F=%f;*');
    vars.eta = sscanf(fgetl(fileID), 'eta=%f;*');
    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.S = sscanf(fgetl(fileID), 'S=%f;*');
    vars.Us = sscanf(fgetl(fileID), 'Us=%f;*');
    vars.xi = sscanf(fgetl(fileID), 'xi=%f;*');
    vars.Gamma = sscanf(fgetl(fileID), 'Gamma=%f;*');
    vars.q_rad = sscanf(fgetl(fileID), 'q_rad=%f;*');
    
    tline = fgetl(fileID);          % blank line
    for i = 1:5
        tline = fgetl(fileID);      % header line
    end

    vars.uS_shape = sscanf(fgetl(fileID), 'u_shape=%[a-z];');
    vars.H_uS_const = sscanf(fgetl(fileID), 'H_uS_const=%f;');
    vars.a_qbo = sscanf(fgetl(fileID), 'a_qbo=%f;');
    vars.b_qbo = sscanf(fgetl(fileID), 'b_qbo=%f;');
    vars.c_qbo = sscanf(fgetl(fileID), 'c_qbo=%f;'); 

    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.z_cirrus = sscanf(fgetl(fileID), 'z_cirrus=%f;');
    vars.upsilon = sscanf(fgetl(fileID), 'upsilon=%f;');
    vars.u_advect = sscanf(fgetl(fileID), 'u_advect=%f;');

    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.H = sscanf(fgetl(fileID), 'H=%f;'); 
    vars.H_s = sscanf(fgetl(fileID), 'H_s=%f;');
    vars.p_s = sscanf(fgetl(fileID), 'p_s=%f;');
    vars.p_t = sscanf(fgetl(fileID), 'p_t=%f;');
    vars.T_bl = sscanf(fgetl(fileID), 'T_bl=%f;');
    vars.cp_d = sscanf(fgetl(fileID), 'cp_d=%f;');
    vars.Rd = sscanf(fgetl(fileID), 'Rd=%f;');
    vars.g = sscanf(fgetl(fileID), 'g=%f;');
    vars.Lv = sscanf(fgetl(fileID), 'Lv=%f;');
    vars.a = sscanf(fgetl(fileID), 'a=%f;');
    vars.beta = sscanf(fgetl(fileID), 'beta=%f;');
    
    tline = fgetl(fileID);          % blank line
    for i = 1:4
        tline = fgetl(fileID);      % header line
    end
    init_from_file = sscanf(fgetl(fileID), 'init_from_file=%[a-z];');
    vars.init_from_file = strcmpi(init_from_file,'true');
    vars.init_file = sscanf(fgetl(fileID), 'init_file=%[A-Za-z/.0-9_:];');

    tline = fgetl(fileID);          % blank line
    for i = 1:4
        tline = fgetl(fileID);      % header line
    end
    vars.Nx = sscanf(fgetl(fileID), 'Nx=%f;');
    vars.Ny = sscanf(fgetl(fileID), 'Ny=%f;');
    vars.Np = sscanf(fgetl(fileID), 'Np=%f;');
    vars.Nz_lower = sscanf(fgetl(fileID), 'Nz_lower=%f;');
    vars.Nz_upper = sscanf(fgetl(fileID), 'Nz_upper=%f;');
    vars.res_edge = sscanf(fgetl(fileID), 'res_edge=%f;');
    vars.z_top = sscanf(fgetl(fileID), 'z_top=%f;');
    vars.y_sponge = sscanf(fgetl(fileID), 'y_sponge=%f;');
    vars.z_sponge = sscanf(fgetl(fileID), 'z_sponge=%f;');
    vars.r_sponge = sscanf(fgetl(fileID), 'r_sponge=%f;');    

    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.T = sscanf(fgetl(fileID), 'T=%f;');
    vars.n_steps = sscanf(fgetl(fileID), 'n_steps=%f;');

    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    plot_steps = sscanf(fgetl(fileID), 'plot_steps=%[a-z];');
    vars.plot_steps = strcmpi(plot_steps,'true');
    vars.rescale_threshold = sscanf(fgetl(fileID), 'rescale_threshold=%f;');
    vars.save_state_step = sscanf(fgetl(fileID), 'save_state_step=%f;');
    vars.base_dir = sscanf(fgetl(fileID), 'base_dir=%[A-Za-z/.0-9_:];');
    vars.exp_name = sscanf(fgetl(fileID), 'exp_name=%[A-Za-z/.0-9_:];');

    tline = fgetl(fileID);          % blank line
    tline = fgetl(fileID);          % header line
    vars.T_analyze_start = sscanf(fgetl(fileID), 'T_analyze_start=%f;');
    vars.T_analyze_end = sscanf(fgetl(fileID), 'T_analyze_end=%f;');
end
