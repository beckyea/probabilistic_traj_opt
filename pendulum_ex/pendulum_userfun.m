function [F] = pendulum_userfun(X)
%PENDULUM_USERFUN Defines nonlinear part 
global N Q R;
x = X(1:N*2);
u = X(N*2+1:(N*3-1));
h = X(end);
F = zeros(N*2+1, 1);
F(1) = pendulum_lw(x, u, h, Q, R);
for i = 1:N-1
    x_i = [x(2*i-1) x(2*i)];
    x_ip1 = [x(2*i+1) x(2*i+2)];
    u_i = u(i);
    F(1) = F(1) + cost(x_i, u_i, h);
    [theta_fun, thetadot_fun] = f_h(x, u_i, h, x_ip1);
    F(2*i) = theta_fun;
    F(2*i+1) = thetadot_fun;
end

end

function [c] = cost(~, ~, h)
c = h;
end

function [theta_fun, thetadot_fun] = f_h(x, u, h, x_ip1)
global g L m;
    theta_i = x(1);
    thetadot_i = x(2);
    theta_ip1 = x_ip1(1);
    thetadot_ip1 = x_ip1(2);
    theta_fun = theta_ip1 - (theta_i + thetadot_i*h);
    thetadot_fun = thetadot_ip1 - (thetadot_i - g*sin(theta_i)/L + u*h/(m*L*L));
end
