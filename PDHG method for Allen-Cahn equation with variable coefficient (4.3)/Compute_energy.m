% Compute the free energy given the numerical solution u based on formula (4.5)


function [energy] = Compute_energy(L, N_x, mu, a, b, u)

   
    h_x = L / N_x;

    nodes_rshift = (0.5 + (0 : (N_x-1))) * h_x;    
    sinx_rshift = sin(nodes_rshift);
    siny_ushift = transpose(sin(nodes_rshift));
    
    sigma = (1 + mu / 2 * (sinx_rshift.*sinx_rshift + siny_ushift.*siny_ushift));

    u_rshift = u;
    u_rshift(:, 2:N_x) = u(:, 1:N_x-1);
    u_rshift(:, 1) = u(:, N_x);

    u_ushift(2:N_x, :) = u(1:N_x-1, :);
    u_ushift(1, :) = u_ushift(N_x, :);

    nabla_u_sq = (u_rshift - u).*(u_rshift - u) + (u_ushift - u).*(u_ushift - u);

    W_u = (u.*u - 1).*(u.*u - 1) / 4; 
    energy = sum( a / 2 * sigma .* nabla_u_sq + b * W_u * h_x * h_x, "All");

end


