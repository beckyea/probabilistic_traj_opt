function x = LQR_x_d(t)
global thetas thetadots dt
i_lower = floor(t/dt)+1;
t_lower = (i_lower-1)*dt;
global N
if t/dt < N - 1
    theta = (thetas(i_lower+1) - thetas(i_lower)) * (t-t_lower)/dt + thetas(i_lower);
    thetadot = (thetadots(i_lower+1) - thetadots(i_lower)) * (t-t_lower)/dt + thetadots(i_lower);
else 
    theta = thetas(i_lower);
    thetadot = thetadots(i_lower);
end
x = [theta thetadot]';
end