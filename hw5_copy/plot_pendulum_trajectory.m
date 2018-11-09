function plot_pendulum_trajectory(t, x)

figure(1)
th = x(:,1)';

B = 2;
xrange = [min(th) - B, max(th) + B];
tic

p = params();
L = p.l;
pend = .1;
pennred = [149,0,26]/256;
px = L*sin(th);
py = - L*cos(th);


stale = .01;
tic

i = 1;

while i<=numel(t)
    start = toc;
    hold off;
     
    plot(4*xrange,[0 0], 'k', 'LineWidth',3)
    hold on;
    plot([0, L*sin(th(i))],[0, -L*cos(th(i))], 'k', 'LineWidth',2);
    
    rectangle('Position',[L*sin(th(i))-pend/2,-L*cos(th(i))-pend/2,pend,pend],...
        'Curvature',[1,1], 'FaceColor',pennred,'EdgeColor','k','LineWidth',2);
    plot(px(1:i), py(1:i), 'g','LineWidth',2);
    axis equal;
    %xlim([x(i) - 2, x(i) + 2]);
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    %legend('x_d(t)','x(t)');
    xlabel('y');
    ylabel('z');
    titl = sprintf('Cart-Pole Trajectory, $t =  %.2f $',t(i));
    title(titl,'Interpreter','latex');
    
    
    compu = toc - start;
    stale_i = max(stale,compu*2);
    next_i = find(t >= start + stale_i);
    if numel(next_i) < 1
        if i < numel(t)
            i = numel(t);
        else
            break;
        end
    else
        i = next_i(1);
    end
    pause(t(i) - toc);

end