function xdot = compute_cart_pole_dyn(x, t, u)
%function to compute the continuous time dynamics of the cartpole for use
%with ode45
global Mp Mc g L mu_k Nc_sign 
theta =  x(1);
thetadot = x(2);
y = x(3);
ydot = x(4);

%friction force Not distinguising between static and kinetic friction
%right now
% Ff = mu_k*(Mc + Mp)*g*sign(ydot); %is this the correct normal force
% 
% M = [         Mp*L^2, Mp*L*cos(theta);
%      Mp*L*cos(theta),          Mc+Mp];
% C = [           Mp*g*L*sin(theta); 
%      -Mp*L*sin(theta)*thetadot^2];
% B = [0; u - Ff];
% 
% qddot = M\(B-C);
% thetaddot = qddot(1);
% yddot = qddot(2);

% new try with florian equations
% thetaddot_num = g*sin(theta) + cos(theta)*((-u - Mp*L*thetadot^2*(sin(theta) ...
%                 + mu_k*sign(Nc_sign*ydot)*cos(theta)))/(Mc + Mp) + mu_k*g*sign(Nc_sign*ydot));
% thetaddot_den = L*(4/3 - (Mp*cos(theta))/(Mc + Mp)*(cos(theta) - mu_k*sign(Nc_sign*ydot)));
% thetaddot = thetaddot_num/thetaddot_den;
S = Nc_sign*sign(ydot);
thetaddot_num = -(Mp+Mc)*g*sin(theta)-u*cos(theta)-Mp*L*thetadot^2*sin(theta) + ...
    mu_k*S*((Mp+Mc)*cos(theta)+Mp*L*thetadot^2*cos(theta));
thetaddot_den = L*(Mp+Mc-Mp*cos(theta)^2+Mp*sin(theta)*cos(theta)*mu_k*S);
thetaddot = thetaddot_num/thetaddot_den;

Nc = (Mc+Mp)*g + Mp*L*(thetaddot*sin(theta) + thetadot^2*cos(theta));

if sign(Nc) ~= Nc_sign
    %recompute thetaddot
    Nc_sign = sign(Nc);
    S = Nc_sign*sign(ydot);
    thetaddot_num = -(Mp+Mc)*g*sin(theta)-u*cos(theta)-Mp*L*thetadot^2*sin(theta) + ...
        mu_k*S*((Mp+Mc)*cos(theta)+Mp*L*thetadot^2*cos(theta));
    thetaddot_den = L*(Mp+Mc-Mp*cos(theta)^2+Mp*sin(theta)*cos(theta)*mu_k*S);
    thetaddot = thetaddot_num/thetaddot_den;
end
Nc = (Mc+Mp)*g + Mp*L*(thetaddot*sin(theta)+thetadot^2*cos(theta));

yddot = (u-mu_k*Nc*S+Mp*L*(thetadot^2*sin(theta)-thetaddot*cos(theta)))/(Mc+Mp);

xdot = [thetadot; thetaddot; ydot; yddot];
