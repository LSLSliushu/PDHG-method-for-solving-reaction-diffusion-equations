% This function applies PDHG method to solve the time-implicit schemes of the Allen-Cahn equation (4.3) with a mobility term on time interval [0, T_total]. 
% We devide [0, T_total] into M intervals, and apply PDHG method to solve for numerical solution on each interval which consists of N_t subintervals.
%
% Try:
% [UT_total, computation_time, accumul_computn_time, Energy_decay] = PDHG_PDE_solver(0.01, 2*pi, 256, 5.0, 500, 1, 1000, 0.95, 0.55, 1, 100);


function [UT_total, computation_time, accumul_comptn_time, Energy_decay] = PDHG_PDE_solver(epsilon_0, L, N_x, T_total, M, N_t, Iter_number, tau_P, tau_U, omega, plotperiod)
% UT_total: numerical solution on [0, T_total], N_x by N_x by M*N_t array
% computation_time: computation time (s)
% accumul_computn_time: accumulated computation time (s), M by 1 array

      
      T = T_total / M;

      mu = 5.0;
      a = epsilon_0;
      b = 1 / epsilon_0;

      hx = L / N_x;
      x_node = transpose(hx * (0 : N_x-1));
      y_node = hx * (0 : N_x-1);
      U_initial = 0.5 * (cos(4 * x_node) + cos(4 * y_node)); 

      UT_total = zeros(N_x, N_x, M * N_t);
      computation_time = 0;
      accumul_comptn_time = zeros(M, 1);
      Energy_decay = zeros(M * N_t, 1);
      
      for num = 1 : M

          fprintf("Computing on the time subinterval [%f, %f]\n", (num-1) * T, num * T);   
   
          %%%%%%%%%%%%%% Apply PDHG algorithm to compute the numerical solution on each time interval %%%%%%%%%%%%%%
          [UT, time] = PDHG_PDE_solver_on_interval_with_initial( epsilon_0, U_initial, L, N_x, T, N_t, Iter_number, tau_P, tau_U, omega );
          UT_total(:, :, (num-1)*N_t+1 : num*N_t) = UT;
          U_initial = UT(:, :, N_t);

          computation_time = computation_time + time;
          accumul_comptn_time(num) = computation_time;

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
                filename = sprintf('[Allen-Cahn with mobility term equation (4.3)] heatmap t=%f .fig', num * T);
                savefig(gcf, filename);
                close();

          end

      end

      %%%%%%%%%%%%%%%% compute and plot the free energy vs time %%%%%%%%%%%%%%
      for t = 1 : M * N_t
          Energy_decay(t) = Compute_energy(L, N_x, mu, a, b, UT_total(:, :, t));
      end      

      ht = T / N_t;
      t_axis = (1 : M * N_t) * ht;
      figure
      plot(t_axis, Energy_decay, 'b');
      title("Free energy decay plot (Implicit scheme solved via PDHG)");
      xlabel('time');
      ylabel('Free energy');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)] Free energy decay.fig'); 
      savefig(gcf, filename);
      close();

      figure
      plot(log(t_axis), log(Energy_decay), 'r');
      title("Free energy decay log-log plot (Implicit scheme solved via PDHG)");
      xlabel('log(time)');
      ylabel('log(Free energy)');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)] log-log plot of free energy decay.fig'); 
      savefig(gcf, filename);
      close();

end



