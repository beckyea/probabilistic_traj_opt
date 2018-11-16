clear
h = 1;
theta1 = 0;
u1 = 2; u2 = 2.5; u3 = 1.5;
thetadot1 = .5;
m = 1;
L = 1;
g = 9.81;
N = 4;

theta2 = theta1 + thetadot1 * h;
thetadot2 = thetadot1 + u1*h/(m*L*L) - g*h*sin(theta1)/L;

theta3 = theta2 + thetadot2 * h;
thetadot3 = thetadot2 + u2*h/(m*L*L) - g*h*sin(theta2)/L;

theta4 = theta3 + thetadot3 * h;
thetadot4 = thetadot3 + u3*h/(m*L*L) - g*h*sin(theta3)/L;

X = [theta1; thetadot1; theta2; thetadot2; theta3; thetadot3;...
 theta4; thetadot4; u1; u2; u3; h]

pendulum_userfun(X)

thetas = zeros(1,N);
t = 0:h:h*(N-1);
for i = 1:N
    thetas(i) = X(2*i-1);
end
pendulum_plot(t,thetas)