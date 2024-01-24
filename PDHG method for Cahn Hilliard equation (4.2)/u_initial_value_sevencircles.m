% COmpute the initial value of equation (4.2)


function [u_initial] = u_initial_value_sevencircles(epsilon, L, n_x)
   N_x = n_x;
   h_x = L/n_x;   

   x1 = pi/2;
   y1 = pi/2;
   r1 = pi/5;

   x2 = pi/4;
   y2 = 3*pi/4;
   r2 = 2*pi/15;
   
   x3 = pi/2;
   y3 = 5*pi/4;
   r3 = 2*pi/15;

   x4 = pi;
   y4 = pi/4;
   r4 = pi/10;

   x5 = 3*pi/2;
   y5 = pi/4;
   r5 = pi/10;

   x6 = pi;
   y6 = pi;
   r6 = pi/4;

   x7 = 3*pi/2;
   y7 = 3*pi/2;
   r7 = pi/4;

   x_node = transpose(h_x * (0:N_x-1));
   y_node = h_x * (0:N_x-1);
   u1 = sqrt((x_node - x1).*(x_node - x1) + (y_node - y1).*(y_node - y1)) - r1;
   m1 = CH_mollifier_function(u1, epsilon);

   u2 = sqrt((x_node - x2).*(x_node - x2) + (y_node - y2).*(y_node - y2)) - r2;
   m2 = CH_mollifier_function(u2, epsilon);

   u3 = sqrt((x_node - x3).*(x_node - x3) + (y_node - y3).*(y_node - y3)) - r3;
   m3 = CH_mollifier_function(u3, epsilon);

   u4 = sqrt((x_node - x4).*(x_node - x4) + (y_node - y4).*(y_node - y4)) - r4;
   m4 = CH_mollifier_function(u4, epsilon);

   u5 = sqrt((x_node - x5).*(x_node - x5) + (y_node - y5).*(y_node - y5)) - r5;
   m5 = CH_mollifier_function(u5, epsilon);

   u6 = sqrt((x_node - x6).*(x_node - x6) + (y_node - y6).*(y_node - y6)) - r6;
   m6 = CH_mollifier_function(u6, epsilon);

   u7 = sqrt((x_node - x7).*(x_node - x7) + (y_node - y7).*(y_node - y7)) - r7;
   m7 = CH_mollifier_function(u7, epsilon);

   u_initial = m1 + m2 + m3 + m4 + m5 + m6 + m7 - 1;

   
end

