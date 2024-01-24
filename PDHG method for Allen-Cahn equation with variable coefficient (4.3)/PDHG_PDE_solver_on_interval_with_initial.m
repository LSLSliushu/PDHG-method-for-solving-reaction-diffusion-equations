% Use PDHG algorithm to solve the implicit scheme on a time interval with
% length equals to T consisting of N_t subintervals.
% 
% As mentioned in section 2.3 of the manuscript, G-prox PDHG method (i.e., (2.17) with the U update replaced by (2.18)) is equivalent to the proposed preconditioned PDHG method (2.15).
% We are implementing the G-prox PDHG algorithm here, since it is more convenient in coding.


function [UT, computn_time] = PDHG_PDE_solver_on_interval_with_initial(epsilon_0, U_initial_2D, L, N_x, T, N_t, Iter_number, tau_P, tau_U, omega)
   
    tic    

    fprintf("========================  \n");
    fprintf("The PDHG algorithm starts.\n");

    a = epsilon_0;
    b = 1 / epsilon_0;
    mu = 5.0;
    ave_sgm = 1 + mu * 1/2;
    ht = T/N_t;
    
    
    epsilon = 0.1;  
    tol = 1e-6;
 
    
    UT = zeros(N_x, N_x, N_t);
    PT = zeros(N_x, N_x, N_t);
    for j = 1 : N_t
       UT(:, :, j) = U_initial_2D;
    end


    %%%%%%%%%%%%% PDHG method starts here %%%%%%%%%%%%%%
    for n = 1 : Iter_number

        if mod(n, 50) == 0
            disp(n);
        end

        %%%%%%%%%%%%%% update the dual variable P %%%%%%%%%%%%%%
        PT_old = PT;
        F_UT = FFT_FUT(UT, U_initial_2D, L, T, N_t, N_x, a, b, mu);
        G_inv_F_UT = inv_G_AllenCahn_with_mobility(F_UT, T, N_t, L, N_x, a, b, ave_sgm);
        PT = 1 / (1 + tau_P * epsilon) * (PT + tau_P * G_inv_F_UT);
        
        %%%%%%%%%%%%%% extrapolation step %%%%%%%%%%%%%
        tilde_PT = PT + omega * (PT - PT_old);

        %%%%%%%%%%%%%% update the primal variable U %%%%%%%%%%%%%%%%
        DF_transpose_PT = FFT_DF_transpose_P(UT, tilde_PT, ht, N_t, L, N_x, a, b, mu);
        UT = UT - tau_U * DF_transpose_PT;
        
        
        %%%%%%%%%%%%%% terminate the iteration if || Res ||_\infty < tol %%%%%%%%%%%%%%
        res_UT = norm(F_UT/ht, 'fro');
        max_res_UT = max(abs(F_UT/ht), [], "all");
        if (max_res_UT < tol)
            break;
        end

        if mod(n-1, 50) == 0
           fprintf('iteration %d.  \n', n);
           fprintf('L2 residue of UT=');
           disp(res_UT);
           fprintf('max residual of UT =');
           disp(max_res_UT);
        end          

    end

    fprintf("The PDHG method finishes after %d iterations\n", n);
    fprintf("========================  \n");


    computn_time = toc;


end







