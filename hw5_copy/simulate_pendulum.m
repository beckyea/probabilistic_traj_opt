function simulate_pendulum
% SIMULATE_CARTPOLE simulates a trajectory for the cartpole system
N = 60;
dt = 6/N;
nx = 2;
nu = 1;

x_0 = [0; 0];
x_f = [pi; 0];

z = find_swingup_trajectory(x_0, x_f, N, dt);

u_t = 0:dt:(dt*(N-1));
u = zeros(nu,0);
x = zeros(nx,0);
for i=1:N
   x_i_inds = (1:nx) + (nx + nu) * (i - 1);
   u_i_inds = (1:nu) + nx * i + nu * (i - 1);
   
   x(:,i) = z(x_i_inds);
   u(:,i) = z(u_i_inds);
end

[K, S] = lqr(A_mat(dt),B_mat(dt),eye(nx),eye(nu));
        
opts = odeset('MaxStep', 0.1,'RelTol',1e-4,'AbsTol',1e-4);

[t, x] = ode45(@(t,x) dynamics(t, x, u, dt, K, S), [0 dt*(N)*1.5], x_0, opts);

plot_pendulum_trajectory(t, x);

t = u_t;

figure(2)
hold on
plot(t, x)
plot(t, u_t)
xlabel('Time')
ylabel('Value')
legend('\theta', '\theta_{dot}', 'u')
title('States and Inputs versus Time')
hold off

end

function dx = dynamics(t, x, u, dt, K, S)
    x_star = [pi; 0];
    if (x-x_star)'*S*(x-x_star) < 3
        u_x = - K*(x-x_star);
    else
        i = ceil(t/dt);
        if i == 0
            i = 1;
        end
        if i >= numel(u)
            i = numel(u) - 1;
        end
        lambda = mod(t,dt)/dt;
        % time is between u(i) and u(i+1)
        u_x = u(i)*(1-lambda) + u(i+1)*lambda;
    end
    
    dx = f(x, u_x);
end
