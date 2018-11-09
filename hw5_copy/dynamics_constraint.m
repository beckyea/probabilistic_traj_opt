function h_i = dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt)
%DYNAMICS_CONSTRAINT(xm, um, xp, up, dt) computes the vector contstraint
%   h_i = \dot{s}_i(\Delta t/2) - f(s_i(\Delta t/2), 0.5(u_{i+1}+u{i})) = 0
%
%   @param x_i: x_i, the state at the start of the interval.
%   @param u_i: u_i, the input at the start of the interval.
%   @param x_ip1: x_{i+1}, the state at the end of the interval.
%   @param u_ip1: u_{i+1}, the input at the end of the interval.
%   @param dt: \Delta t, the duration of the interval
%
%   @output h_i: quantity from above expresion that should be 0


f_i = f(x_i, u_i);
f_ip1 = f(x_ip1, u_ip1);

s_i_c = (x_ip1+x_i)/2 - dt*(f_ip1-f_i)/8;
sdot_i_c = 3*(x_ip1-x_i)/(2*dt)-(f_ip1+f_i)/4;

h_i = sdot_i_c - f(s_i_c, (u_ip1+u_i)/2);

end