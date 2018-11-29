function [A,B] = cart_pole_cont_dyn_lin() 
%% compute continous dynamics linearizations
syms theta thetadot theta y ydot y h w u g Mp Mc L theta_D thetadot_D y_D ydot_D u_D mu_k
thetaddot = (-(Mp+Mc)*g*sin(theta)-cos(theta)*(u+Mp*L*thetadot^2*sin(theta)))/...
        (L*(Mp+Mc)-Mp*L*cos(theta)*cos(theta));
Nc = 0; S = 0;
yddot = (u-mu_k*Nc*S+Mp*L*(thetadot^2*sin(theta)-thetaddot*cos(theta)))/(Mc+Mp);
xdot = [thetadot; thetaddot; ydot; yddot];
x = [theta; thetadot; y; ydot];
x_D = [theta_D; thetadot_D; y_D; ydot_D];
e = x - x_D;
v = u - u_D;

A = simplify(jacobian(xdot, x));
B = simplify(jacobian(xdot, u));

end