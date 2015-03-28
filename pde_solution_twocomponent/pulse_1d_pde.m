%% Will McFadden (wmcfadden)
function [c, f, s] = pulse_1d_pde(x,t,u,dudx,Da,Dr,l,L,m0,K,n,kon_a,koff_a,kon_r,koff_r)
    %% initialization
    a = u(1);
    r = u(2);
    v = u(3);
    dadx = dudx(1);
    drdx = dudx(2);
    dvdx = dudx(3);
    
    m = m0*(a.^n)./(K^n+a.^n);
        
    c = [1;1;0];
    
    f = [Da*dadx - a*v; ...
         Dr*drdx - r*v; ...
         l^2*dvdx + m];
    
    s = [kon_a-koff_a*a; ...
         kon_r-koff_r*r; ...
         -v];
    
end