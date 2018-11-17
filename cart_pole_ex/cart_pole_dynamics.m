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
