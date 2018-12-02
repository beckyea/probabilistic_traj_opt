
function u = LQR_u_d(t)
global N us dt
i_lower = floor(t/dt)+1;
t_lower = (i_lower-1)*dt;
if t/dt < N - 2
    u = (us(i_lower+1) - us(i_lower)) * (t-t_lower)/dt + us(i_lower);
elseif t/dt < N - 1
    u = us(i_lower);
else
    u = 0;
end
end