%% Will McFadden (wmcfadden)
function [c, f, s] = active_1d_pde(x,t,u,dudx,D,l,L,m0,K,n,kon,koff)
    %% initialization
    a = u(1);
    v = u(2);
    dadx = dudx(1);
    dvdx = dudx(2);
    
    m = m0*(a.^n)./(K^n+a.^n);
        
    c = [1;0];
    
    f = [D*dadx - a*v; ...
        l^2*dvdx + m];
    
    s = [kon-koff*a; ...
        -v];
    
end