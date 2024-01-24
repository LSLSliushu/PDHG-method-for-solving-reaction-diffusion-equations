function [vec] = laplace_vector(L, n_x) 
   h_x = L / n_x;
   N_x = n_x;
   vec = zeros(N_x, 1);
   vec(1) = -2 / (h_x * h_x);
   vec(2) = 1 / (h_x * h_x);
   vec(N_x) = 1 / (h_x * h_x);

end