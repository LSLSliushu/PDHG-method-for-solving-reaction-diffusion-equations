% Use FFT to compute X^{-1}U.
% Here:
% X = I - ht (a * Lap - b * c * I);
% U is N_x * N_x.


function [inv_X_U] = inv_X(ht, L, N_x, U, a, b)

     c = 2;
 
     Lap_vec = laplace_vector(L, N_x);
     FFT_Lap_vec = fft(Lap_vec);
     FFT_Lap_x_tensor = reshape(FFT_Lap_vec, [N_x, 1, 1]);
     FFT_Lap_y_tensor = reshape(FFT_Lap_vec, [1, N_x, 1]);

     FFT_X = 1 - ht * ( a * (FFT_Lap_x_tensor + FFT_Lap_y_tensor) - b * c);
     
     FFT_U = fft2(U);
     FFT_inv_X_U = FFT_U./FFT_X;
     inv_X_U = ifft2(FFT_inv_X_U);

     inv_X_U = real(inv_X_U);

end