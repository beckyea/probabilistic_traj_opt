function f = f(x,u)

% Defined Values
p = params();

theta = x(1,:);
thetadot = x(2,:);

f = [thetadot;
     - p.g/p.l*sin(theta)+u/(p.m*p.l^2)];
% pdot = x(3,:);
% thetadot = x(4,:);
% theta = x(2,:);
% t2 = cos(theta);
% t3 = sin(theta);
% t4 = t2.^2;
% t5 = t4-2.0;
% t6 = 1.0./t5;
% t7 = thetadot.^2;
% f = [pdot;thetadot;t6.*(u.*1.0e2+t2.*t3.*9.81e2+t3.*t7.*1.0e2).*(-1.0./1.0e2);t6.*(t3.*9.81e2+t2.*u.*5.0e1+t2.*t3.*t7.*5.0e1).*(1.0./5.0e1)];
