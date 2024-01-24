% This function applies PDHG method with adaptive interval sizes to solve the time-implicit schemes of the Allen-Cahn equation (4.3)  on time interval [0, T_total]. 
% We start computing on the first interval with length equals t1. We devide this interval into N_t subintervals (each having length equal to t1 / N_t), 
% we apply PDHG method to compute on these N_t subintervals. Then, we move on to compute on the second interval with length t2.
% We set t2 = t1/2 if the previous PDHG takes more than MAX_PDHG_ITER_NUM steps to converge (i.e., ||Res|| < tol), 
% and t2 = 1.1 * t1 if the PDHG takes no more than MAX_PDHG_ITER_NUM steps to converge. 
% We repeat such process until the time reaches T_total.


% Try:
% [UT_total, computation_time, accumul_computn_time, t_list, accumul_time, iter_num_list] = PDHG_PDE_solver_adaptive_time_stepsize(0.01, 2*pi, 256, 20.0, 0.01, 1, 4000, 0.95, 0.55, 1, 0.08, 400, 50);


function [U_T_total, computation_time, accumul_comptn_time, t_list, accumul_t_list, iter_num_list] = PDHG_PDE_solver_adaptive_time_stepsize(epsilon_0, L, N_x, T_total, t1, N_t, Iter_number, tau_P, tau_U, omega, MAX_TIME_INTERVAL_LENGTH, MAX_PDHG_ITER_NUM, plotperiod)
% U_T_total: numerical solution at (around) T_total, N_x by N_x array
% computation_time: computation time (s)
% accumul_computn_time: accumulated computation time (s), M by 1 array
% t_list: a list of adptive size of each time interval 
% accumul_t_list: accumulated size of time intervals
% iter_num_list: iterations needed for the PDHG method on each time interval

      hx = L / N_x;
      x_node = transpose(hx * (0 : N_x-1));
      y_node = hx * (0 : N_x-1);
      U_initial = 0.5 * (cos(4 * x_node) + cos(4 * y_node)); 
      
      computation_time = 0;

      MAX_SIZE = 2 * floor(T_total / t1);
      accumul_comptn_time = zeros(MAX_SIZE, 1);  % preallocate the array
      t_list = zeros(MAX_SIZE, 1);  % preallocate the array
      accumul_t_list = zeros(MAX_SIZE, 1);  % preallocate the array
      iter_num_list = zeros(MAX_SIZE, 1);  % preallocate the array
      
      time_interval_length = t1;

      num = 0;
      while true

          num = num+1;
          t_list(num) = time_interval_length;
          if num > 1
              accumul_t_list(num) = accumul_t_list(num-1) + time_interval_length;
          else
              accumul_t_list(num) = time_interval_length;
          end

          fprintf("Computing on the time subinterval [%f, %f]\n", accumul_t_list(num), accumul_t_list(num)+time_interval_length);   
   
          %%%%%%%%%%%%%% Apply PDHG algorithm to compute the numerical solution on each time interval %%%%%%%%%%%%%%
          [U, time, PDHG_iter_num] = PDHG_PDE_solver_on_interval_with_initial( epsilon_0, U_initial, L, N_x, time_interval_length, N_t, Iter_number, tau_P, tau_U, omega );
          U_initial = U;

          computation_time = computation_time + time;
          accumul_comptn_time(num) = computation_time;

          iter_num_list(num) = PDHG_iter_num;

          %%%%%%%%%%%%%% adjust the length of time interval %%%%%%%%%%%
          if PDHG_iter_num <= MAX_PDHG_ITER_NUM
              if 1.1 * time_interval_length < MAX_TIME_INTERVAL_LENGTH
                  time_interval_length = 1.1 * time_interval_length;
              end
          else
              time_interval_length = 0.5 * time_interval_length;
          end
          
          %%%%%%%%%%%%%%%% plot the heat map of the numerical solution every plotperiod intervals %%%%%%%%%%%%%%%
          if mod(num, plotperiod) == 0

                UT_2D = U;
                figure
                
                X = (0:N_x-1)*hx;
                Y = (0:N_x-1)*hx;
                imagesc(X, Y, UT_2D);
                xlim([0, L]);
                ylim([0, L]);
                colorbar;
                set(gca, 'YDir', 'normal');
                filename = sprintf('[Allen-Cahn with mobility term equation (4.3)] heatmap t=%f .fig', accumul_t_list(num));
                savefig(gcf, filename);
                close();

          end

          if accumul_t_list(num) >= T_total
              break;
          end

      end

      M = num;

      U_T_total = U;
      %%%%%%%%%%%%%%%% plot the heat map of the numerical solution at (around) T_total %%%%%%%%%%%%%%%
      UT_2D = U_T_total;
      figure
      X = (0:N_x-1)*hx;
      Y = (0:N_x-1)*hx;
      imagesc(X, Y, UT_2D);
      xlim([0, L]);
      ylim([0, L]);
      colorbar;
      set(gca, 'YDir', 'normal');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)] heatmap t=%f .fig', accumul_t_list(M));
      savefig(gcf, filename);
      close();


      %%%%%%%%%%%%%%%% plot adaptive time interval lengths vs accumulated time %%%%%%%%%%%%%%%%%%%%%%%      
      figure
      used_t_list = t_list(1:M);
      used_accumul_t_list = accumul_t_list(1:M);
      plot(used_accumul_t_list, used_t_list);        
      xlabel('accumulated time  (physical time of equation) ');
      ylim([0 0.1]);
      ylabel('adaptive time interval length');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)]  adaptive time interval length - accumulated time   plot,  T_total = %f.fig', T_total);
      savefig(gcf, filename);
      close();
      

      
      %%%%%%%%%%%%%%%% plot PDHG iteration number on each time interval vs accumulated time %%%%%%%%%%%%
      figure
      used_iter_num_list = iter_num_list(1:M);
      plot(used_accumul_t_list, used_iter_num_list);        
      xlabel('accumulated time  (physical time of equation) ');
      ylabel('PDHG iteration number');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)]  PDHG iteration number on each interval - accumulated time   plot,  T_total = %f.fig', T_total);
      savefig(gcf, filename);
      close();



      %%%%%%%%%%%%%%%% plot computation time (CPU time) vs accumulated time (physical time of equation) %%%%%%%%%%%%%%%%%%%%%%%
      figure
      used_accumul_computn_time = accumul_comptn_time(1:M);
      plot( used_accumul_t_list, used_accumul_computn_time );
      xlabel(' accumulated time  (physical time of equation) ');
      ylabel('accumulated computation time  (CPU time, in seconds)');
      filename = sprintf('[Allen-Cahn with mobility term equation (4.3)]  accumulated computation time (CPU time) - accumulated time (physical time)   plot,  T_total = %f.fig', T_total);
      savefig(gcf, filename);
      close();





end



