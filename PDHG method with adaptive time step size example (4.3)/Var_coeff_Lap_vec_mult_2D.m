% UT is N_x by N_x by N_t
% Use central difference scheme to approximate ∇∙(σ∇U).


function [L_UT] = Var_coeff_Lap_vec_mult_2D(L, N_t,N_x, UT, a)

    h_x = L / N_x;
    
    nodes = (0 : (N_x - 1)) * h_x;
    nodes_rshift = (0.5 + (0 : (N_x-1))) * h_x;
    nodes_lshift = (-0.5 + (0 :(N_x-1))) * h_x;

    sinx = sin(nodes);
    siny = transpose(sin(nodes));
    sinx_rshift = sin(nodes_rshift);
    siny_ushift = transpose(sin(nodes_rshift));
    sinx_lshift = sin(nodes_lshift);
    siny_dshift = transpose(sin(nodes_lshift));
    
    U_rightshift = zeros(N_x, N_x, N_t);
    U_rightshift(1:N_x-1, :, :) = UT(2:N_x, :, :);
    U_rightshift(N_x, :, :) = UT(1, :, :);

    U_leftshift = zeros(N_x, N_x, N_t);
    U_leftshift(2:N_x, :, :) = UT(1:N_x-1, :, :);
    U_leftshift(1,:,:) = UT(N_x, :, :);

    U_upshift = zeros(N_x, N_x, N_t);
    U_upshift(:, 1:N_x-1, :) = UT(:, 2:N_x, :);
    U_upshift(:, N_x, :) = UT(:, 1, :);

    U_downshift = zeros(N_x, N_x, N_t);
    U_downshift(:, 2:N_x, :) = UT(:, 1:N_x-1, :);
    U_downshift(:, 1, :) = UT(:, N_x, :);

    L_xx_UT = ((1 + a / 2 * (sinx_rshift.*sinx_rshift + siny.*siny)).*(U_rightshift - UT) - (1 + a / 2 * (sinx_lshift.*sinx_lshift + siny.*siny)).*(UT - U_leftshift)) / (h_x * h_x);
    L_yy_UT = ((1 + a / 2 * (sinx.*sinx + siny_ushift.*siny_ushift)).*(U_upshift - UT) - (1+a/2*(sinx.*sinx + siny_dshift.*siny_dshift)).*(UT - U_downshift)) / (h_x * h_x);
    L_UT = L_xx_UT + L_yy_UT;
    
end
    