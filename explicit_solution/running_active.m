% Will McFadden (wmcfadden)
% plots velocities for active fluid with various fluid parameters and
% concentration profiles

s = 100;
cc = jet(s);


q = [50 1 0.5 100];
x = (1:0.1:100)';
m = x./x;
m0 = x./x-x./x;
v0=0;
vd=0;


for i=1:100
    
    m = m0+i*exp(-((x-max(x)/2).^2))/100;
        
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    hold on
    plot(x, m,'color',cc(floor(0.25*s),:),'LineWidth',3)
    ylim([-1 1])
    drawnow
end

for i=1:49
    
    q(1) = q(1)-1;    
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    hold on
    plot(x, m,'color',cc(floor(0.25*s),:),'LineWidth',3)
    ylim([-1 1])
    drawnow
end


for i=1:50
    
    m = m0+exp(-((x-max(x)/2).^2)/i);
        
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    hold on
    plot(x, m,'color',cc(floor(0.25*s),:),'LineWidth',3)
    ylim([-1 1])
    drawnow
end

for i=1:100
    
    q(1) = q(1)+0.1;    
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    hold on
    plot(x, m,'color',cc(floor(0.25*s),:),'LineWidth',3)
    ylim([-1 1])
    drawnow
end
