% Use FFT to compute X^{-1}U.
% Here:
% X = I_N_x*N_x - ht △(ε_0^2 △ - cI + ε_0^2I)(ε_0^2 △ - cI),
% here △ indicates discretization of Laplace operator w.r.t. periodic b.c.
% c = W''(±1) = 2.
%
% U is N_x * N_x.


function [inv_X_U] = inv_X(ht, L, N_x, U, epsilon_0)

     c  =  2;
 
     Lap_vec = laplace_vector(L, N_x);
     FFT_Lap_vec = fft(Lap_vec);
     FFT_Lap_x_tensor = reshape(FFT_Lap_vec, [N_x, 1, 1]);
     FFT_Lap_y_tensor = reshape(FFT_Lap_vec, [1, N_x, 1]);

     FFT_X = 1 - ht * (FFT_Lap_x_tensor + FFT_Lap_y_tensor) .* ( epsilon_0 * epsilon_0 * (FFT_Lap_x_tensor + FFT_Lap_y_tensor) -  c  + epsilon_0 * epsilon_0 ) .* (epsilon_0 * epsilon_0 * (FFT_Lap_x_tensor + FFT_Lap_y_tensor) - c);
     
     FFT_U = fft2(U);
     FFT_inv_X_U = FFT_U./FFT_X;
     inv_X_U = ifft2(FFT_inv_X_U);

     inv_X_U = real(inv_X_U);

end