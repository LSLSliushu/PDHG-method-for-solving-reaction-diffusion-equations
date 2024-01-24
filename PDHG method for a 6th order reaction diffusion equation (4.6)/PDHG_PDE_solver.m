% This function applies PDHG method to solve the time-implicit schemes of the 6th order reaction diffusion equation (4.6) on time interval [0, T_total]. 
% We devide [0, T_total] into M intervals, and apply PDHG method to solve for numerical solution on each interval which consists of N_t subintervals.
%
% Try:
% [UT_total, computation_time, accumul_computn_time] = PDHG_PDE_solver(0.18, 2*pi, 256, 10.0, 1000, 1, 1500, 0.95, 0.5, 1, 50);


function [UT_total, computation_time, accumul_computational_time] = PDHG_PDE_solver(epsilon_0, L, N_x, T_total, M, N_t, Iter_number, tau_P, tau_U, omega, plotperiod)
% UT_total: numerical solution on [0, T_total], N_x by N_x by M*N_t array
% computation_time: computation time (s)
% accumul_computn_time: accumulated computation time (s), M by 1 array

      T = T_total / M;

      UT_total = zeros(N_x, N_x, M * N_t);
      computation_time = 0;
      accumul_computational_time = zeros(N_t, 1);

      % set up the initial value
      hx = L / N_x;
      x_node = transpose(hx * (0:N_x-1));
      y_node = hx * (0:N_x-1);
      U_0 = 2 * exp( sin(x_node) + sin(y_node) - 2 ) + 2.2 * exp( - sin(x_node) - sin(y_node) - 2 ) - 1;
      U_0 = transpose(U_0);


      U_initial = U_0;

      for num = 1 : M

          fprintf("Computing on the time subinterval [%f, %f]\n", (num-1) * T, num * T); 
   
          %%%%%%%%%%%%%% Apply PDHG algorithm to compute the numerical solution on each time interval %%%%%%%%%%%%%%
          [UT, time] = PDHG_PDE_solver_on_interval_with_initial(epsilon_0, U_initial, L, N_x, T, N_t, Iter_number, tau_P, tau_U, omega);
          UT_total(:, :, (num-1)*N_t+1 : num*N_t) = UT;
          U_initial = UT(:, :, N_t);

          computation_time = computation_time + time;
          accumul_computational_time(num) = computation_time;

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
                filename = sprintf('[6th order reaction diffusion equation (4.6)] heatmap t=%f .fig', num * T);
                savefig(gcf, filename);
                close();

          end

      end

      fprintf("Computation time (s): ");
      disp(computation_time);



end



