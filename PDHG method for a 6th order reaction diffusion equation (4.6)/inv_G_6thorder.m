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
% X = I_N_x*N_x - ht △(ε_0^2 △ - 2I + ε_0^2I)(ε_0^2 △ - 2I), here △ indicates discretization of the Laplace operator w.r.t. periodic b.c..
% X is N_x^2 by N_x^2.


function [inv_G_Y] = inv_G_6thorder(Y, T, N_t, L, N_x, epsilon_0)

     ht = T / N_t;

     inv_M_Y = zeros(N_x, N_x, N_t);
     inv_M_Y(:, :, 1) = inv_X(ht, L, N_x, Y(:, :, 1), epsilon_0);
     for i = 2:N_t
         inv_M_Y(:, :, i) = inv_X(ht, L, N_x, inv_M_Y(:, :, i-1) + Y(:, :, i), epsilon_0); 
     end
     
     inv_G_Y = zeros(N_x, N_x, N_t);
     inv_G_Y(:, :, N_t) = inv_X(ht, L, N_x, inv_M_Y(:, :, N_t), epsilon_0);
     for i = (N_t-1):-1:1
         inv_G_Y(:,:, i) = inv_X(ht, L, N_x, inv_G_Y(:, :, i+1) + inv_M_Y(:, :, i), epsilon_0);
     end

end
