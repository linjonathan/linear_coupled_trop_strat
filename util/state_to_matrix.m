function [state_matrix] = state_to_matrix(uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y, uS_yz, vS_yz, phiS_yz, eta_y, qT_y)
    state_trop = cat(1, uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y);
    state_matrix = cat(1, state_trop, uS_yz, vS_yz, phiS_yz, eta_y, qT_y);
end