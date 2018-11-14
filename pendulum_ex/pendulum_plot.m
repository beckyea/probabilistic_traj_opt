function pendulum_plot(t, th)

figure(1)
clf

B = 2;
xrange = [min(th) - B, max(th) + B];
tic

global L;
pend = .1;
px = L*sin(th);
py = - L*cos(th);


stale = .01;
tic

i = 1;

while i<=numel(t)
    start = toc;
    hold off;
     
    plot(4*xrange,[0 0], 'k', 'LineWidth',2)
    hold on;
    plot([0, L*sin(th(i))],[0, -L*cos(th(i))], 'k', 'LineWidth',1);
    
    rectangle('Position',[L*sin(th(i))-pend/2,-L*cos(th(i))-pend/2,pend,pend],...
        'Curvature',[1,1], 'FaceColor',[0, 0.7, 0.8],'EdgeColor','k','LineWidth',1);
    plot(px(1:i), py(1:i), 'g','LineWidth',1);
    axis equal;
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    xlabel('y');
    ylabel('z');
    titl = sprintf('Pendulum Trajectory, $t =  %.2f $',t(i));
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
    pause(t(i) - toc+ 1);

end