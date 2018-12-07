clc
syms m g L theta_i thetadot_i u_i theta_ip1 thetadot_ip1 h real


theta_fun = theta_ip1 - (theta_i + thetadot_i*h);
thetadot_fun = thetadot_ip1 - (thetadot_i - g*h*sin(theta_i)/L + u_i*h/(m*L*L));

dyn_i = [theta_fun; thetadot_fun];

x_derivs = simplify([diff(dyn_i, theta_i),   diff(dyn_i, thetadot_i),...
                     diff(dyn_i, theta_ip1), diff(dyn_i, thetadot_ip1)]);
                     
u_derivs = [diff(dyn_i, u_i)];

h_derivs = [diff(dyn_i, h)];

%% X Derivs
A_x = [ -1,  0, 1, 0;
         0, -1, 0, 1];

G_x = [ 0,                  -h, 0, 0;
       (g*h*cos(theta_i))/L, 0, 0, 0];

A_u = [ 0; 0];
G_u = [ 0; -h/(L^2*m)];

G_h = [-thetadot_i; (g*sin(theta_i))/L - u_i/(L^2*m)];
    

simplify(x_derivs-A_x-G_x)
simplify(u_derivs-A_u-G_u)
simplify(h_derivs-G_h)