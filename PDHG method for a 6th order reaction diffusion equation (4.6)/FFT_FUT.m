% Computes the function F(UT) defined in (2.4).
% UT is N_x * N_x * N_t.


function [F_UT] = FFT_FUT(UT, U_initial_2D, L, T, N_t, N_x,  epsilon_0)

    ht = T / N_t;
    rhs = sixth_order_Lap_vec_mult_2D(L, N_x, UT, epsilon_0);
    translate_forward_UT = zeros(N_x, N_x, N_t);
    if N_t > 1
       translate_forward_UT(:, :, 2:N_t) = UT(:, :, 1:N_t-1);
    end    
    translate_forward_UT(:, :, 1) = U_initial_2D;

    F_UT = UT - translate_forward_UT - ht * rhs;

end



