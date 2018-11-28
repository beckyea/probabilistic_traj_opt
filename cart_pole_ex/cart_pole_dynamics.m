%cart pole dynamics script
syms g Mp Mc theta thetadot y ydot L u
%% calculate accelerations
M = [Mc + Mp , Mp*L*cos(theta); 
     Mp*L*cos(theta) , Mp*L^2];
 
C = [-Mp*L*sin(theta)*thetadot^2; Mp*g*L*sin(theta)];

B = [1; 0];

q = [y theta]';
qdot = [ydot; thetadot];

qddot = M\(B*u - C); %NOTE: outputs in the form [yddot; thetaddot]

f = [thetadot; qddot(2); ydot; qddot(1)]; % NOTE: outputs in the form [thetadot; thetaddot; ydot; yddot]
%% calculate the discrete dynmaics
syms theta_i thetadot_i theta_ip1 y_i ydot_i y_ip1 h
thetaddot_i = subs(qddot(2), theta, theta_i); 
thetaddot_i = subs(thetaddot_i, thetadot, thetadot_i);
yddot_i = subs(qddot(1), theta, theta_i);
yddot_i = subs(yddot_i, thetadot, thetadot_i);

x_i = [theta_i; thetadot_i; y_i; ydot_i];
xdot_i = [thetadot_i; thetaddot_i; ydot_i; yddot_i];
x_ip1 = x_i + xdot_i * h; %discrete dynamics

%% calculate discrete linearizations
A = simplify(jacobian(x_ip1, x_i));
B = simplify(jacobian(x_ip1, u));

%% try to use the other dynamics with w instead of modeled friction
syms theta_i thetadot_i theta_ip1 y_i ydot_i y_ip1 h w u_i g Mp Mc L
% thetaddot_i = (g*sin(theta_i) + cos(theta_i)*((-u -Mp*L*(sin(theta_i) + w))/(Mc + Mp))) / ...
%     (L*(4/3 - Mp*cos(theta_i)/(Mc + Mp)*(cos(theta_i) - w)));
% Nc = (Mc + Mp)*g - Mp*L*(thetaddot_i*sin(theta_i) + thetadot_i^2*cos(theta_i));
% yddot_i = (u + Mp*L*(thetadot_i^2*sin(theta_i) - thetaddot_i*cos(theta_i)) - w*Nc) / ...
%     (Mc + Mp);
thetaddot_i = (-(Mp+Mc)*g*sin(theta_i)-cos(theta_i)*(u_i-w+Mp*L*thetadot_i^2*sin(theta_i)))...
    /(L*(Mp+Mc)-Mp*L*cos(theta_i)*cos(theta_i));

yddot_i = (u_i - w + Mp*L*(thetadot_i^2*sin(theta_i)-thetaddot_i*cos(theta_i)))/(Mc+Mp);

x_i = [theta_i; thetadot_i; y_i; ydot_i];
xdot_i =[thetadot_i; thetaddot_i; ydot_i; yddot_i];
x_ip1 = x_i + xdot_i * h;
A = simplify(jacobian(x_ip1, x_i));
B = simplify(jacobian(x_ip1, u_i));
G = simplify(jacobian(x_ip1, w));
A = simplify(subs(A, w, 0));
B = simplify(subs(B, w, 0));
G = simplify(subs(G, w, 0));

