function [uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y, uS_yz, vS_yz, phiS_yz, eta_y, qT_y] = matrix_to_state(state_matrix, vars)
    uBT_y = state_matrix(1, :);
    uBC_y = state_matrix(2, :);
    vBT_y = state_matrix(3, :);
    vBC_y = state_matrix(4, :); 
    sT_y = state_matrix(5, :);
    smT_y = state_matrix(6, :);
    idx_offset = 7; idx_var = vars.Nz;
    uS_yz = state_matrix(idx_offset:(idx_offset+idx_var-1), :);
    vS_yz = state_matrix((idx_offset+idx_var):(idx_offset+idx_var*2-1), :);
    phiS_yz = state_matrix((idx_offset+idx_var*2):(idx_offset+idx_var*3-1), :);
    eta_y = state_matrix(end-1, :);
    qT_y = state_matrix(end, :);
end