function [z, Aeq, beq, lb, ub, z0] = find_swingup_trajectory(x_0, x_f, N, dt)
%FIND_SWINGUP_TRAJECTORY(x_0, x_f, N, dt) executes a direct collocation
%optimization program to find an input sequence to drive the cartpole
%system from x_0 to x_f.
%
%   @param x_0: the state at the start of the trajectory
%   @param x_f: the state at the emd of the trajectory
%   @param N: number of state/input samples
%   @param dt: \Delta t, the duration of the interval between samples
%
%   @output z: decision variable vector containing the x_i and u_i
%   @output Aeq: matrix from linear constrant Aeq z = beq
%   @output beq: vector from linear constrant Aeq z = beq
%   @output lb: lower bound from box constraint lb <= z <= ub
%   @output ub: upper bound from box constraint lb <= z <= ub
%   @output z0: initial guess for z

  nx = 2;
  nu = 1;
  
  % TODO: Add constraints to Aeq, beq to enforce starting at x_0 and ending
  % at x_f
  x_0_inds = 1:nx;
  x_f_inds = x_0_inds + (N - 1) * (nx + nu);
  Aeq = zeros(2*nx, N * (nx + nu));
  beq = zeros(2*nx, 1);
  
  for i = 1:nx
      x_0_ind = x_0_inds(i);
      x_f_ind = x_f_inds(i);
      Aeq(i, x_0_ind) = 1;
      Aeq(i+nx, x_f_ind) = 1;
      beq(i, 1) = x_0(i);
      beq(i+nx, 1) = x_f(i);
  end
  
  M = 50;
  
  % TODO: Add bounding box constraints u_1 \in [-M,M], u_2 \in [-M,M]
  lb = -inf(N * (nx + nu),1);
  ub = inf(N * (nx + nu),1);
  for i=1:N
      u_i_inds = (1:nu) + nx * i + nu * (i - 1);
      lb(u_i_inds) = -M;
      ub(u_i_inds) = M;
  end
  
  % Make initial guess for z
  z0 = zeros(N * (nx + nu), 1);
  for i=1:N
      x_i_inds = (1:nx) + (nx + nu) * (i - 1);
      z0(x_i_inds,1) = (x_f-x_0)*(i-1)/(N-1) + x_0;
  end
  
  options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Display','iter');
  problem.objective = @(z) trajectory_cost(z, N, nx, nu, dt);
  
  
  problem.x0 = z0;
  problem.options = options;
  problem.nonlcon = @(z) all_constraints(z, N, nx, nu, dt);
  problem.solver = 'fmincon';
  problem.Aeq = Aeq;
  problem.beq = beq;
  problem.lb = lb;
  problem.ub = ub;

  z = fmincon(problem);
end