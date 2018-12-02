x_e0 = 0;

global Q_N thetas thetadots us dt N 
pendulum_globals();
thetas = dirtrel_thetas; 
thetadots = dirtrel_thetadots; 
us = dirtrel_u;
dt = h;
tf = dt * (N-1);

L0sq = chol(Q_N)';
mask = tril(true(2,2));
L0 = reshape(L0sq,[4 1]);

opts = odeset('MaxStep',0.01);
[t, Lmat] = ode45(@dLdtminus, [0 tf], L0, opts);
t = tf - t;
t = flip(t);
Lmat = flipud(Lmat);
Lmat = spline(t,Lmat');

[t, x] = ode45(@(t,x) f(x, feedback(t, x, ppval(Lmat, t))), [0 tf], LQR_x_d(0) + x_e0);

% get Us from ODE call
% u = zeros(size(t,1)-1,1);
% for i = 1:size(t,1)-1
%     u(i) = feedback(t,x,Lmat);
% end
%     

pendulum_plot(t, x(:,1), [1; zeros(size(x-2,1),1)])
%pendulum_plot(0:dt:(N-1)*dt, thetas, dirtrel_u)


function Ldotminus = dLdtminus(t,L)
    Lsq = reshape(L,[2 2]);
    dLsq = LQR_Ldot(t, Lsq);
    Ldotminus = -reshape(dLsq,[4 1]);
end

function u = feedback(t, x, Lmat)
    global R
    pendulum_globals();
    B = [0 1/(m*L*L)]';
    Lsq = reshape(Lmat,[2 2]);
    S = Lsq*Lsq';
    u = LQR_u_d(t);
    u_e = -inv(R)*B'*S*(x-LQR_x_d(t));
    u = u + u_e;
end

function xdot = f(xs, us)
theta = xs(1); 
thetadot = xs(2);
u = us(1);
global m g L
xdot = [thetadot; -g/L*sin(theta)+u/(m*L*L)];
end