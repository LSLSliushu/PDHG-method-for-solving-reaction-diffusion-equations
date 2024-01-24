% Compute the front position (distance from the zero level set of the numerical solution to the center (L/2, L/2) of the square [0, L]^2) of the numerical solution U.
% U is N by N


function [x] = front_position(U, N, L)
   h_x = L / N;
   x = 0;
   
   half_N = floor(N/2)+1;
   for k = half_N : N-1
       if U(half_N, k) * U(half_N, k+1) <= 0
           x = (k-1) * h_x - L/2 + U(half_N, k) / (U(half_N, k)-U(half_N, k+1)) * h_x;
           break;
       end
   end

end



