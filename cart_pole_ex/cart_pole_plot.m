function cart_pole_plot(t, th, y, u)

set(gcf,'color','white')

B = 2;
xrange = [min(th) - B, max(th) + B];
tic

global L N;
h = .2;
w = .4;
pend = .1;
py = y + L*sin(th);
pz = - L*cos(th);

blue = [0, 0.7, 0.8];
purp = [0.4, 0.2, 0.4];
gray = [0.5, 0.5, 0.5];


stale = .01;
tic

i = 1;

while i<=numel(t)
    start = toc;
    subplot(2,1,1)
    hold off;
     
    plot(4*xrange,[-h/2 -h/2], 'k', 'LineWidth',2)
    hold on;
    rectangle('Position',[y(i)-w/2, -h/2, w, h],'FaceColor',purp,'EdgeColor',purp,...
    'LineWidth',3)

    plot([y(i), y(i) + L*sin(th(i))],[0, -L*cos(th(i))], 'k', 'LineWidth',1);
    
    rectangle('Position',[y(i)+L*sin(th(i))-pend/2,-L*cos(th(i))-pend/2,pend,pend],...
        'Curvature',[1,1], 'FaceColor',blue,'EdgeColor', blue,'LineWidth',2);
    plot(py(1:i), pz(1:i), 'Color',gray,'LineWidth',1);
    axis equal;
    xlim([-2.5 2.5]);
    ylim([-1.5 1.5]);
    xlabel('y');
    ylabel('z');
    titl = sprintf('Trajectory, t =  %.2f',t(i));
    title(titl);
    
        subplot(2,1,2)
    if i< N
        plot(t(1:i), u(1:i), 'Color', blue, 'LineWidth', 1.25)
    else
        plot(t, [u; 0], 'Color', blue, 'LineWidth', 1.25)
    end
    title('Control');
    ylabel('u')
    xlabel('t')
    ylim([-max(-min(u),max(u)) max(-min(u),max(u))]);
    xlim([min(t) max(t)]);
    
    
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