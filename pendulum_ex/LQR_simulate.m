x_e0 = [0 0];

global thetas thetadots us dt N Q_N_test

% Define LQR Controller
Q = [10 0; 0 10];
R = 0.001;

pendulum_globals();
thetas = dirtrel_thetas; 
thetadots = dirtrel_thetadots; 
us = dirtrel_u;
dt = h;
tf = dt * (N-1);

L0sq = chol(Q_N_test)';
mask = tril(true(2,2));
L0 = reshape(L0sq,[4 1]);

opts = odeset('MaxStep',0.01);
[t, Lmat] = ode45(@dLdtminus, [0 tf], L0, opts);
t = tf - t;
t = flip(t);
Lmat = flipud(Lmat);
Lmat = spline(t,Lmat');
[t, x] = ode45(@(t,x) f(x, feedback(t, x, ppval(Lmat, t))), [0 tf], LQR_x_d(0), opts);

% get Us from ODE call
u = zeros(size(t,1)-1,2);
for i = 1:size(t,1)-1
    u(i,:) = feedback(t(i),x(i,:),ppval(Lmat, t(i)));
end
    

pendulum_plot(t, x(:,1), u(:,1), 0);
%pendulum_plot(0:dt:(N-1)*dt, thetas, dirtrel_u)


function Ldotminus = dLdtminus(t,L)
    Lsq = reshape(L,[2 2]);
    dLsq = LQR_Ldot(t, Lsq);
    Ldotminus = -reshape(dLsq,[4 1]);
end

function u = feedback(t, x, Lmat)
    global R_test m_test
    pendulum_globals();
    B = [0 1/(m_test*L*L)]';
    Lsq = reshape(Lmat,[2 2]);
    S = Lsq*Lsq';
    u = LQR_u_d(t);
%     x_d = LQR_x_d(t)';
%     x_e = (x-LQR_x_d(t))';

    u_e = -inv(R_test)*B'*S*(x-LQR_x_d(t));

%     disp("x: [" + x(1) + " " + x(2) + "]"  + " x_d: [" + x_d(1) + " " + x_d(2) + "]" + " x_e: [" + x_e(1) + " " + x_e(2) + "]");
% 
%     disp("u: " + u + " u_e: " + u_e + " new_u:" + (u+u_e));
    u = u + u_e;
end

function xdot = f(xs, us)
theta = xs(1); 
thetadot = xs(2);
u = us(1);
global m_test g L
xdot = [thetadot; -g/L*sin(theta)+u/(m_test*L*L)];
% disp("theta: " + theta + " thetadot: " + thetadot + " u: " + u)
% disp("thetad: " + xdot(1) + " thetadotd: " + xdot(2))

end