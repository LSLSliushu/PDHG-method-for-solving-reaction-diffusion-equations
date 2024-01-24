

function [vec] = laplace_vector(L, N_x)

   h_x = L / N_x;
   vec = zeros(N_x, 1);
   vec(1) = -2 / (h_x * h_x);
   vec(2) = 1 / (h_x * h_x);
   vec(N_x) = 1 / (h_x * h_x);

end