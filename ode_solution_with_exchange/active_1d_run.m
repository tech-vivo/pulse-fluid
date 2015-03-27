% Will McFadden (wmcfadden)
% simulates the 1d active fluid model in the presence of binding and unbinding


T = 1000;       %total time to simulate 
samps = 100;     %number of timepoints to sample solution

D = 0.5;        %diffusion coefficient
L = 100;        %domain size
l = 5;          %ratio of viscosity to friction
c0 = 1;         %"equilibrium concentration"
koff = 0.001;   %off rate for reaction
kon = c0*koff;  %on rate (must be derived from off rate and equilibrium conc)

%these parameters define how to convert from concentration to active stress
m0 = 10;        
K = c0;
n = 1;



%initialize discrete positions and corresponding concentration
x = linspace(0,L,200)';
c = c0*(ones(size(x))-0.01*(rand(size(x))-0.5)-0.01*cos(2*pi*x/L));

%integrate
[t, c] = ode23tb(@active_1d_ode,linspace(0,T,samps),c,odeset('NonNegative',1:length(c),'RelTol',1e-3),x,m0,K,n,D,l,L,kon,koff);
x = x';

%plot the solutions
for i = 1:size(c,1)
    plot3(t(i)*ones(size(x)),x,c(i,:));
    hold on
end

%%make a movie plotting the concentration profile through time
% figure;
% movind = 1;
% for i = 1:size(c,1)
%     m = m0*(c(i,:).^n)./(K^n+c(i,:).^n);
%     Gr = (cosh((L+x(1)-x)/l)-cosh((x(1)-x)/l))/2/l^2/(cosh(L/l)-1);
%     v(1) = trapz(x,Gr.*m);
%     Gl = (cosh((L-x(end)+x)/l)-cosh((x(end)-x)/l))/2/l^2/(cosh(L/l)-1);
%     v(length(x)) = -trapz(x,Gl.*m);
%     for ind = 2:length(x)-1
%         Gr = (cosh((L+x(ind)-x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
%         Gl = (cosh((L-x(ind)+x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
%         v(ind) = trapz(x(ind:end),Gr(ind:end).*m(ind:end)) -  trapz(x(1:ind),Gl(1:ind).*m(1:ind));
%     end
%     
%     plot(x,c(i,:));
%     ylim([0 10*c0]);
%     
%     drawnow
%     mov(movind,:) = getframe;
%     movind = movind+1;
% end
