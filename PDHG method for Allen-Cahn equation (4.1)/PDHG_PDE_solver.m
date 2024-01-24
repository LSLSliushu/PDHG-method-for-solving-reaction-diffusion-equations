% This function applies PDHG method to solve the time-implicit schemes of the Allen-Cahn equation (4.1) on time interval [0, T_total]. 
% We devide [0, T_total] into M intervals, and apply PDHG method to solve for numerical solution on each interval which consists of N_t subintervals.
%
% Try:
% [UT_total, computation_time, accumul_computn_time, front_position_list ] = PDHG_PDE_solver(0.01, 0.5, 256, 2.5, 2500, 1, 1000, 0.95, 0.55, 1, 50);


function [UT_total, computation_time, accumul_computn_time, front_position_list] = PDHG_PDE_solver( epsilon_0, L, N_x, T_total, M, N_t, Iter_number, tau_P, tau_U, omega, plotperiod)
% UT_total: numerical solution on [0, T_total], N_x by N_x by M*N_t array
% computation_time: computation time (s)
% accumul_computn_time: accumulated computation time (s), M by 1 array
% front_position_list: recording the front positions, M*N_t by 1 array

    T = T_total / M;
    hx = L / N_x;
    r = 0.2;
    mid_L = L / 2;

    x_node = transpose(hx * (0:N_x-1));
    y_node = hx * (0:N_x-1);
    U_initial = 2.0 * ((x_node-mid_L).*(x_node-mid_L) + (y_node-mid_L).*(y_node-mid_L) <= r*r) - 1;

    UT_total = zeros(N_x, N_x, M * N_t);      
    front_position_list = zeros(M*N_t, 1);
    
    computation_time = 0;
    accumul_computn_time = zeros(M, 1);

    for num = 1 : M

        fprintf("Computing on the time subinterval [%f, %f]\n", (num-1) * T, num * T);
   
        %%%%%%%%%%%%%% Apply PDHG algorithm to compute the numerical solution on each time interval %%%%%%%%%%%%%%
        [UT, time, front_position_list_on_interval] = PDHG_PDE_solver_on_interval_with_initial(epsilon_0, U_initial, L, N_x, T, N_t, Iter_number, tau_P, tau_U, omega);
        UT_total(:, :, (num-1)*N_t+1 : num*N_t) = UT;
        U_initial = UT(:, :, N_t);

        front_position_list((num-1)*N_t+1:num*N_t) = front_position_list_on_interval;
        
        computation_time = computation_time + time;
        accumul_computn_time(num) = computation_time;

        %%%%%%%%%%%%%%%% plot the heat map of the numerical solution every plotperiod intervals %%%%%%%%%%%%%%
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
            filename = sprintf('[Allen-Cahn equation (4.1)] heatmap t=%f .fig',  num * T);
            savefig(gcf, filename);
            close();
    
        end

    end

    fprintf("Computation time (s): ");
    disp(computation_time);

    %%%%%%%%%%%%%%%% plot the front positions vs time %%%%%%%%%%%%%%
    ht = T_total / (M * N_t);
    figure
    t_node = ht: ht : ht * M * N_t;
    plot(t_node, front_position_list, 'b-');
    xlabel("Time");
    ylabel("Front position");
    filename = sprintf('[Allen-Cahn equation (4.1)] front position (L = %f, N_x = %d, epsilon_0=%f).fig', L, N_x, epsilon_0);
    savefig(gcf, filename);
    close();

end



