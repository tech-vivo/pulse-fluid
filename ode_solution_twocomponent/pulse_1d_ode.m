%% Will McFadden (wmcfadden)

function db = pulse_1d_ode(t,b,x,m0,K,n,Da,Dr,l,L,kon_a,koff_a,kon_r,koff_r)
    %% initialization
    t
    a = b(1:length(b)/2);
    r = b(length(b)/2+1:end);
    
    dx = x(2)-x(1);             % space discretization size
    da = zeros(size(a));        % initialize dc/dt array to return
    dr = zeros(size(r));        % initialize dr/dt array to return
    v = zeros(size(x));         % initial velocity array
    
    m = m0*(a.^n)./(K^n+a.^n);  % convert concentration into active stress
    
    
    %% Compute the velocity from active stress m with fluid length scale l 
    
    Gr = (cosh((L+x(1)-x)/l)-cosh((x(1)-x)/l))/2/l^2/(cosh(L/l)-1);        
    v(1) = trapz(x,Gr.*m);
    Gl = (cosh((L-x(end)+x)/l)-cosh((x(end)-x)/l))/2/l^2/(cosh(L/l)-1);
    v(length(x)) = -trapz(x,Gl.*m);
    for ind = 2:length(x)-1
        Gr = (cosh((L+x(ind)-x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
        Gl = (cosh((L-x(ind)+x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
        v(ind) = trapz(x(ind:end),Gr(ind:end).*m(ind:end)) -  trapz(x(1:ind),Gl(1:ind).*m(1:ind));
    end
    
    %%  use velocity and diffusion coefficient to compute changes in concentrations
    
    for ind = 1:length(v)
        left = ind-1;
        if(left<1)
            left = length(v);
        end
        right = ind+1;
        if(right>length(v))
            right = 1;
        end
        jL = (a(left)+a(ind))*(v(left)+v(ind))/4 - Da*(a(ind)-a(left))/dx;
        jR = (a(right)+a(ind))*(v(right)+v(ind))/4 - Da*(a(right)-a(ind))/dx;
        da(ind) = da(ind) - (jR - jL)/dx;
        da(left) = da(left) - jL/dx;
        da(right) = da(right) + jR/dx;
        
        jL = (r(left)+r(ind))*(v(left)+v(ind))/4 - Dr*(r(ind)-r(left))/dx;
        jR = (r(right)+r(ind))*(v(right)+v(ind))/4 - Dr*(r(right)-r(ind))/dx;
        dr(ind) = da(ind) - (jR - jL)/dx;
        dr(left) = da(left) - jL/dx;
        dr(right) = da(right) + jR/dx;
    end
    
    %% add the change in concentration due to simple reaction kinetics
    da = da + (kon_a - koff_a*a);
    dr = dr + (kon_r - koff_r*r);
    db = [da; dr];
end