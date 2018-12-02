function [dL] = LQR_Ldot(t, Lmat)
x = LQR_x_d(t);
u = LQR_u_d(t);

theta = x(1);

global m L g Q R

A = [0 1; -g*cos(theta) 0];

B = [0; 1/(m*L*L)];


dL = -.5*Q*inv(Lmat)' - A'*Lmat + 0.5*(Lmat*Lmat')*B*inv(R)*B'*Lmat;

end