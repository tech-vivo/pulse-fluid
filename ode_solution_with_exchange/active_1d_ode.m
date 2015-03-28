%% Will McFadden (wmcfadden)

function dc = active_1d_ode(t,c,x,m0,K,n,D,l,L,kon,koff)
    %% initialization

    dx = x(2)-x(1);             % space discretization size
    dc = zeros(size(c));        % initialize dc/dt array to return
    v = zeros(size(c));         % initial velocity array
    
    m = m0*(c.^n)./(K^n+c.^n);  % convert concentration into active stress
    
    
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
    
    %%  use velocity and diffusion coefficient to compute changes in concentration
    
    for ind = 1:length(v)
        left = ind-1;
        if(left<1)
            left = length(v);
        end
        right = ind+1;
        if(right>length(v))
            right = 1;
        end
        jL = (c(left)+c(ind))*(v(left)+v(ind))/4 - D*(c(ind)-c(left))/dx;
        jR = (c(right)+c(ind))*(v(right)+v(ind))/4 - D*(c(right)-c(ind))/dx;
        dc(ind) = dc(ind) - (jR - jL)/dx;
        dc(left) = dc(left) - jL/dx;
        dc(right) = dc(right) + jR/dx;
    end
    
    %% add the change in concentration due to simple reaction kinetics
    dc = dc + (kon - koff*c);
end