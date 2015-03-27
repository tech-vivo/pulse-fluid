% Will McFadden (wmcfadden)
% plots velocities for passive fluid

s = 100;
cc = jet(s);


q = [1 1 1 1];
x = (1:100)';
m = x./x;
v0=0;
vd=0;
for i=1:100
    v0 = v0+0.01;
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    ylim([0 1])
    drawnow
end


for i=1:25
    
    vd = vd+0.01;
    
    
    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    ylim([0 1])
    drawnow
end

for i=1:150
    q(1) = q(1)+1;

    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    ylim([0 1])
    drawnow
end


for i=1:100
    q(1) = q(1)-1;

    tempdat = {x, m, v0, vd};
    fitdat = {tempdat};
    v = act_flu_fun(q, fitdat);

    plot(x, [v0; v; vd],'color',cc(floor(0.75*s),:),'LineWidth',3);
    ylim([0 1])
    drawnow
end


