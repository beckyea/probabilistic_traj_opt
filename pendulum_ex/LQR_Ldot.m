function [dL] = LQR_Ldot(t, Lmat)
x = LQR_x_d(t);
u = LQR_u_d(t);

theta = x(1);

global m_test L g Q_test R_test

A = [0 1; -g*cos(theta)/L 0];

B = [0; 1/(m_test*L*L)];


dL = -.5*Q_test*inv(Lmat)' - A'*Lmat + 0.5*(Lmat*Lmat')*B*inv(R_test)*B'*Lmat;

end