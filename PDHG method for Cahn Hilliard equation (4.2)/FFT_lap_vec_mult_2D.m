% U is N_x by N_x by N_t


function [Lap_U] = FFT_lap_vec_mult_2D(L, n, U)

   L_vec = laplace_vector(L, n);
   F_L_vec = fft(L_vec);
   F_u_2d = fft2(U);
   M = F_L_vec + transpose(F_L_vec);
   F_Lap_U = M.*F_u_2d;
   Lap_U = ifft2(F_Lap_U);
   Lap_U = real(Lap_U);

end