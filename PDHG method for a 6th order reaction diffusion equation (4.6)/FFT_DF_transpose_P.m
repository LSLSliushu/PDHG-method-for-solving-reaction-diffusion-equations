% Compute DF(UT)^T PT, 
% UT, PT are N_t * N_x * N_x, DF is the Jacobian of F.
% Recall F(UT) = \mathscr{D} UT + ht \mathscr{G} (a \mathscr{L} UT + b f(UT))
% We have DF(UT)^T PT = \mathscr{D}^T PT + ht (Lap PT) .* (-W'''(UT)) .* (- ε_0^2 Lap UT + W'(UT)) + ht (- ε_0^2 Lap + diag(W''(UT))) (ε_0^2 Lap - diag(W''(UT)) + ε_0^2) (Lap PT)
%         \_________/   \______________/      \__________________________________________________/      \________________________________________________________________________/
%           DF_T_PT        DF_T_PT_1                             DF_T_PT_2                                                       DF_T_PT_3     
%
% W'(U) = U^3-U
% W''(U) = 3U^2-1
% W'''(U) = 6U 


function [DF_T_PT] = FFT_DF_transpose_P(UT, PT, ht, N_t, L, N_x, epsilon_0)

    translate_bkwd_PT = zeros(N_x, N_x, N_t);
    translate_bkwd_PT(:, :, 1:N_t-1) = PT(:, :, 2:N_t);
    DF_T_PT_1 = PT - translate_bkwd_PT;
    
    Lap_U = FFT_lap_vec_mult_2D(L, N_x, UT);
    Lap_P = FFT_lap_vec_mult_2D(L, N_x, PT);

    dW_U = UT.*UT.*UT - UT;
    ddW_U = 3 * UT .* UT - 1;
    dddW_U = 6 * UT;

    DF_T_PT_2 = Lap_P .* (-dddW_U) .* (-epsilon_0 * epsilon_0 * Lap_U + dW_U);

    Lap_Lap_P = FFT_lap_vec_mult_2D(L, N_x, Lap_P);
    tilde_Lap_Lap_P = epsilon_0 * epsilon_0 * Lap_Lap_P - ddW_U .* Lap_P + epsilon_0 * epsilon_0 * Lap_P;
    Lap_tilde_Lap_Lap_P = FFT_lap_vec_mult_2D(L, N_x, tilde_Lap_Lap_P);
    DF_T_PT_3 = - epsilon_0 * epsilon_0 * Lap_tilde_Lap_Lap_P + ddW_U .* tilde_Lap_Lap_P;

    DF_T_PT = DF_T_PT_1 + ht * (DF_T_PT_2 + DF_T_PT_3); 

end
    


      