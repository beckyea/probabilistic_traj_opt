clc
syms Mp Mc g L theta_i thetadot_i u_i y_i ydot_i theta_ip1 thetadot_ip1 y_ip1 ydot_ip1 h real

thetaddot_i = (-(Mp+Mc)*g*sin(theta_i)-cos(theta_i)*(u_i+Mp*L*thetadot_i^2*sin(theta_i)))/...
    (L*(Mp+Mc)-Mp*L*cos(theta_i)*cos(theta_i));

yddot_i = (u_i+Mp*L*(thetadot_i^2*sin(theta_i)-thetaddot_i*cos(theta_i)))/(Mc+Mp);

theta_fun    = theta_ip1   -(theta_i   +thetadot_i*h);
thetadot_fun = thetadot_ip1-(thetadot_i+thetaddot_i*h);
y_fun        = y_ip1       -(y_i       +ydot_i*h);
ydot_fun     = ydot_ip1    -(ydot_i    +yddot_i*h);

dyn_i = [theta_fun; thetadot_fun; y_fun; ydot_fun];

x_derivs = simplify([diff(dyn_i, theta_i), diff(dyn_i, thetadot_i), diff(y_i), diff(ydot_i),...
    diff(dyn_i, theta_ip1), diff(dyn_i, thetadot_ip1), diff(y_ip1), diff(ydot_ip1)]);
u_derivs = [diff(dyn_i, u_i)];

h_derivs = [diff(dyn_i, h)]

%% X Derivs

A_x = [ -1, 0, 1, 1, 1, 0, 1, 1;
         0, 0, 0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0];

G_x = [ 0, -h, 0, 0, 0, 0, 0, 0;
       (h*(L*Mp*(2*cos(theta_i)^2-1)*thetadot_i^2-u_i*sin(theta_i)+...
       Mc*g*cos(theta_i)+Mp*g*cos(theta_i)))/(L*(Mc+Mp-Mp*cos(theta_i)^2))-...
       (2*Mp*h*cos(theta_i)*sin(theta_i)*(cos(theta_i)*(L*Mp*sin(theta_i)*thetadot_i^2+u_i)+....
       g*sin(theta_i)*(Mc+Mp)))/(L*(- Mp*cos(theta_i)^2+Mc+Mp)^2),...
       (2*L*Mp*h*thetadot_i*cos(theta_i)*sin(theta_i))/(L*(Mc+Mp)-L*Mp*cos(theta_i)^2)-1, ...
       0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0;
       -(Mp*h*(Mc*g*cos(2*theta_i)-u_i*sin(2*theta_i)-(Mp*g)/2+(Mp*g*cos(2*theta_i))/2+...
       (L*Mp*thetadot_i^2*cos(3*theta_i))/4+L*Mc*thetadot_i^2*cos(theta_i)-...
       (L*Mp*thetadot_i^2*cos(theta_i))/4))/(- Mp*cos(theta_i)^2+Mc+Mp)^2,...
       -(2*L*Mp*h*thetadot_i*sin(theta_i))/(Mp*sin(theta_i)^2+Mc), 0, 0, 0, 0, 0, 0];

A_u = [ 0; 0; 0; 0];
G_u = [ 0; (h*cos(theta_i))/(L*(Mc+Mp)-L*Mp*cos(theta_i)^2);
        0; -(h*((L*Mp*cos(theta_i)^2)/(L*(Mc+Mp)-L*Mp*cos(theta_i)^2)+1))/(Mc+Mp)];

G_h = [ -thetadot_i;
    (cos(theta_i)*(L*Mp*sin(theta_i)*thetadot_i^2+u_i)+g*sin(theta_i)*(Mc+Mp))/(L*(Mc+Mp)-L*Mp*cos(theta_i)^2);
    -ydot_i;
    -(u_i+L*Mp*(thetadot_i^2*sin(theta_i)+(cos(theta_i)*(cos(theta_i)*(L*Mp*sin(theta_i)*thetadot_i^2+u_i)+g*sin(theta_i)*(Mc+Mp)))/(L*(Mc+Mp)-L*Mp*cos(theta_i)^2)))/(Mc+Mp)];
 

simplify(x_derivs-A_x-G_x)
simplify(u_derivs-A_u-G_u)