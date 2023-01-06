function [] = convert_eFlux_to_init(fn_eFlux, fn_init_output)
    eFlux_mat = load(fn_eFlux);
    sig_rigid = eFlux_mat.vars.sig_rigid;
    state_matrix = eFlux_mat.y_structs;
    vars = eFlux_mat.vars;
    save(fn_init_output, 'sig_rigid', 'state_matrix', 'vars');
end