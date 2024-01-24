% This function applies PDHG method to solve the time-implicit schemes of the Cahn-Hilliard equation (4.2) on time interval [0, T_total]. 
% We devide [0, T_total] into M intervals, and apply PDHG method to solve for numerical solution on each interval which consists of N_t subintervals.
%
% Try:
% [UT_total, computation_time, accumul_computn_time] = PDHG_PDE_solver(0.1, 2*pi, 256, 10.0, 1000, 10, 1000, 0.95, 0.55, 1, 100);


function [UT_total, computation_time, accumul_computn_time] = PDHG_PDE_solver(epsilon_0, L, N_x, T_total, M, N_t, Iter_number, tau_P, tau_U, omega, plotperiod)
% UT_total: numerical solution on [0, T_total], N_x by N_x by M*N_t array
% computation_time: computation time (s)
% accumul_computn_time: accumulated computation time (s), M by 1 array

    T = T_total / M;

    hx = L / N_x;
    U_initial = u_initial_value_sevencircles(0.1, L, N_x); 

    UT_total = zeros(N_x, N_x, M * N_t);
    computation_time = 0;
    accumul_computn_time = zeros(M, 1);     

    for num = 1 : M

        fprintf("Computing on the time subinterval [%f, %f]\n", (num-1) * T, num * T);   
        
        %%%%%%%%%%%%%% Apply PDHG algorithm to compute the numerical solution on each time interval %%%%%%%%%%%%%%
        [UT, time] = PDHG_PDE_solver_on_interval_with_initial(epsilon_0, U_initial, L, N_x, T, N_t, Iter_number, tau_P, tau_U, omega);
        UT_total(:, :, (num-1)*N_t+1 : num*N_t) = UT;
        U_initial = UT(:, :, N_t);
        
        computation_time = computation_time + time;
        accumul_computn_time(num) = computation_time;

        %%%%%%%%%%%%%%%% plot the heat map of the numerical solution every plotperiod intervals %%%%%%%%%%%%%%%
        if mod(num, plotperiod) == 0
    
            UT_2D = UT(:, :, N_t);
            figure
               
            X = (0:N_x-1)*hx;
            Y = (0:N_x-1)*hx;
            imagesc(X, Y, UT_2D);
            xlim([0, L]);
            ylim([0, L]);
            colorbar;
            set(gca, 'YDir', 'normal');
            filename = sprintf('[Cahn-Hilliard equation (4.2)] heatmap t=%f .fig',  num * T);
            savefig(gcf, filename);
            close();
    
        end

    end

    fprintf("Computation time (s): ");
    disp(computation_time);


end



