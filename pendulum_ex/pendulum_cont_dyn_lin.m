function [A,B] = pendulum_cont_dyn_lin() 
%% compute continous dynamics linearizations
syms theta thetadot theta u g Mp Mc L theta_D thetadot_D u_D mu_k
thetaddot = -g/Lsin(theta)+u/(m*L*L);

xdot = [thetadot; thetaddot];
x = [theta; thetadot];
x_D = [theta_D; thetadot_D];
e = x - x_D;
v = u - u_D;

A = simplify(jacobian(xdot, x));
B = simplify(jacobian(xdot, u));

end