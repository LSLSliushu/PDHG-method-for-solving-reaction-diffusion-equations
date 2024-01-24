% UT = [U_1, ..., U_T] is N_x by N_x by N_t
% output [\Delta * (\epsilon_0^2 \Delta - (3 * U_1 * U_1 - 1) + \epsilon_0^2) (epsilon_0^2 * \Delta U_1 - ( U_1 * U_1 * U_1 - U_1)), ... ]


function [L_UT] = sixth_order_Lap_vec_mult_2D(L, N_x, UT, epsilon_0)
    
    % f(u) = W'(u) = u^3-u
    f = UT.*UT.*UT - UT;  

    Delta_UT = FFT_lap_vec_mult_2D(L, N_x, UT);
    
    VT = epsilon_0 * epsilon_0 * Delta_UT - f;
    
    Delta_VT = FFT_lap_vec_mult_2D(L, N_x, VT);
    
    % f'(u) = W''(u) = 3u^2-1
    df = 3 * UT .* UT - 1;

    Varcoeff_L_VT = epsilon_0 * epsilon_0 * Delta_VT - df.*VT + epsilon_0 * epsilon_0 * VT;
    
    L_UT = FFT_lap_vec_mult_2D(L, N_x, Varcoeff_L_VT);

end

