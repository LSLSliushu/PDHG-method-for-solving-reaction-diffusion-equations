% Use back substitution to compute G^{-1}Y, Y is N_t * N_x * N_x.
%
% G = MM^T , G^{-1} = M^{-T}M^{-1}. So G^{-1}Y = M^{-T}(M^{-1}Y)
%
% M = [ X            ]
%     [-I X          ]
%     [  -I X        ]
%     [     ... ...  ]
%     [          -I X]
% M is N_tN_x^2 by N_tN_x^2.
%
% X = I_N_x*N_x + ht (a * Lap * Lap - b * c * Lap)
% X is N_x^2 by N_x^2.


function [inv_G_Y] = inv_G_CahnHilliard(Y, T, N_t, L, N_x, a, b)

     ht = T / N_t;

     inv_M_Y = zeros(N_x,N_x,N_t);
     inv_M_Y(:, :, 1) = inv_X(ht, L, N_x, Y(:, :, 1), a, b);
     for i = 2:N_t
         inv_M_Y(:, :, i) = inv_X(ht, L, N_x, inv_M_Y(:, :, i-1) + Y(:, :, i), a, b); 
     end

     inv_G_Y = zeros(N_x, N_x, N_t);
     inv_G_Y(:, :, N_t) = inv_X(ht, L, N_x, inv_M_Y(:, :, N_t), a, b);
     for i = (N_t-1):-1:1
         inv_G_Y(:,:, i) = inv_X(ht, L, N_x, inv_G_Y(:, :, i+1) + inv_M_Y(:, :, i), a, b);
     end

end
