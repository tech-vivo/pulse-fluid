% Will McFadden (wmcfadden)
% draws a pretty picture of dots moving in a flow field

h = 1;
num = 2000;
w = 1;
expand = 400;
v0=0;
vd=0.0;
k=1;
xs = h*rand(num,1);
ys = w*rand(num,1);

lastTime = 20;
distr = exp(-50*(xs-0.5).^2);
xp = [];
yp = [];
for i=1:length(xs)
    if(rand<distr(i))
        xp = [xp xs(i)];
        yp = [yp ys(i)];
    end
end
mov = moviein(lastTime);

xs = xp';
ys = yp';

for i=1:lastTime
    grid = zeros(expand*h,expand*w);
%     grid(1:ceil(max(expand*xs)),:)=0.2;
    for j=1:length(xs)
        x =ceil(expand*xs(j));
        y =ceil(expand*ys(j));
        
        if(x>0&&x<size(grid,1))
            grid(x,y)=1;
        end
    end
    colormap('gray');imagesc(grid');
    mov(:,i)=getframe;
    
    d=w;
    tstep = 0.1;
    v = v0.*exp(xs/k)+(vd-v0*exp(-d/k)).*sinh(xs/k)/sinh(d/k);
    plusterm = [];
    minusterm = [];
    dx = 0.001;
    for k=1:length(xs)
        xh = xs(k):dx:d;
        xl = 0:dx:xs(k);
        
        plusterm = [plusterm dx*sum(cosh((d-xh)/k).*exp(-50*(xh-d/2).^2))];
        minusterm = [minusterm dx*sum(cosh(xl/k).*exp(-50*(xl-d/2).^2))];
    end
    v = v + sinh(xs/k)/sinh(d/k).*plusterm';
    v = v - sinh((d-xs)/k)/sinh(d/k).*minusterm';
    xs = xs +tstep*v;
    if(i==100)
        vd=-vd;
    end
end
grid = zeros(expand*h,expand*w);
   
