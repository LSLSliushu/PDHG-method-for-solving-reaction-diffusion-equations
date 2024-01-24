% Compute DF(UT)^T PT, 
% UT, PT are N_t * N_x * N_x, DF is the Jacobian of F.
% Recall F(UT) = \mathscr{D} UT + ht \mathscr{L} (a \mathscr{L} UT + b f(UT))
% DF(UT) = \mathscr{D} + ht\mathscr{L}(a \mathscr{L} + b diag(f'(UT)))
% Thus, DF(UT)^T PT = \mathscr{D}^T PT + ht (a \mathscr{L} \mathscr{L} + b diag(f'(UT)) \mathscr{L}) PT
%       \_________/   \______________/      \_________________________________________________________/
%         DF_T_PT        DF_T_PT_1                               DF_T_PT_2          
%
% Here:
% \mathscr{L} = [ -Lap                       ]
%               [       -Lap                 ]
%               [             ....           ]
%               [                   -Lap     ]
%               [                        -Lap]  


function [DF_T_PT] = FFT_DF_transpose_P(UT, PT, ht, N_t, L, N_x, a, b)


    translate_bkwd_PT = zeros(N_x, N_x, N_t);
    translate_bkwd_PT(:, :, 1:N_t-1) = PT(:, :, 2:N_t);
    DF_T_PT_1 = PT - translate_bkwd_PT;
    
    Lap_PT = FFT_lap_vec_mult_2D(L, N_x, PT);
    LapLap_PT = FFT_lap_vec_mult_2D(L, N_x, Lap_PT);
    DF_T_PT_2 = a * LapLap_PT - b * (3 * UT.*UT.*Lap_PT - Lap_PT);

    DF_T_PT = DF_T_PT_1 + ht * DF_T_PT_2;


end
    


      