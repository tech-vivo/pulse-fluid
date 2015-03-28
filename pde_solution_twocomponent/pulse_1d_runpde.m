% Will McFadden (wmcfadden)
% simulates the 1d active fluid model in the presence of binding and
% unbinding (using pdepe)

T = 10000;       %total time to simulate 
samps = 200;     %number of timepoints to sample solution
xbins = 200;

Da = 1.2;        %diffusion coefficient
Dr = 1.2;        %diffusion coefficient
L = 100;        %domain size
l = 5;          %ratio of viscosity to friction
a0 = 1;         %"equilibrium concentration"
r0 = 1;         %"equilibrium concentration"
koff_a = 0.008;   %off rate for reaction
kon_a = a0*koff;  %on rate (must be derived from off rate and equilibrium conc)
koff_r = 0.008;   %off rate for reaction
kon_r = r0*koff;  %on rate (must be derived from off rate and equilibrium conc)

%these parameters define how to convert from concentration to active stress
m0 = 5;        
K = a0;
n = 1;

x = linspace(0,L,xbins);
t = linspace(0,T,samps);
%integrate
sol = pdepe(0,@pulse_1d_pde,@ic,@bc,x,t,odeset('Reltol',0.001),Da,Dr,l,L,m0,K,n,kon_a,koff_a,kon_r,koff_r);

%plot the solutions
for i = 1:size(sol,1)
    plot3(t(i)*ones(size(x)),x,sol(i,:,1));
    hold on
end
figure
for i = 1:size(sol,1)
    plot3(t(i)*ones(size(x)),x,sol(i,:,2));
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