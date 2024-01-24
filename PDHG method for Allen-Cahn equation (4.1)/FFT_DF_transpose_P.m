% Compute DF(UT)^T PT, 
% UT, PT are N_t * N_x * N_x, DF is the Jacobian of F.
% Recall F(UT) = \mathscr{D} UT + ht (a \mathscr{L} UT + b f(UT))
% DF(UT) = \mathscr{D} + ht (a \mathscr{L} + b diag(f'(UT)))
% Thus, DF(UT)^T PT = \mathscr{D}^T PT + ht (a \mathscr{L} + b diag(f'(UT))) PT
%       \_________/   \______________/      \_________________________________/
%         DF_T_PT        DF_T_PT_1                       DF_T_PT_2          
% Here 
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
    DF_T_PT_2 = a * (-Lap_PT) + b * (3 * UT.*UT.* PT - PT);

    DF_T_PT = DF_T_PT_1 + ht * DF_T_PT_2;

end
    


      