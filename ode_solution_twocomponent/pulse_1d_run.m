% Will McFadden (wmcfadden)
% simulates the 1d active fluid model in the presence of binding and unbinding
% github hooray

samps = 20;     %number of timepoints to sample solution
xbins = 100;     %spatial discretization size


T = 200;       %total time to simulate 

Da = 0.5;        %diffusion coefficient
Dr = 0.5;        %diffusion coefficient
L = 100;        %domain size
l = 5;          %ratio of viscosity to friction


a0 = 1;         %"equilibrium concentration"
koff_a = 0.001;   %off rate for reaction
kon_a = a0*koff;  %on rate (must be derived from off rate and equilibrium conc)

r0 = 0;         %"equilibrium concentration"
koff_r = 0.001;   %off rate for reaction
kon_r = a0*koff;  %on rate (must be derived from off rate and equilibrium conc)

%these parameters define how to convert from concentration to active stress
m0 = 10;        
K = a0;
n = 1;



%initialize discrete positions and corresponding concentration
x = linspace(0,L,xbins)';
a = a0*(ones(size(x))-0.01*(rand(size(x))-0.5)-0.01*cos(2*pi*x/L));
r = r0*(ones(size(x))-0.01*(rand(size(x))-0.5)-0.01*cos(2*pi*x/L));

b = [a; r];
%integrate
[t, b] = ode23tb(@pulse_1d_ode,linspace(0,T,samps),b,odeset('NonNegative',1:length(b),'RelTol',1e-2),x,m0,K,n,Da,Dr,l,L,kon_a,koff_a,kon_r,koff_r);
a = b(:,1:length(b)/2);
r = b(:,length(b)/2+1:end);
    

%plot the solutions
for i = 1:size(a,1)
    plot3(t(i)*ones(size(x)),x,a(i,:));
    hold on
end
figure
for i = 1:size(r,1)
    plot3(t(i)*ones(size(x)),x,r(i,:));
    hold on
end

%%make a movie plotting the concentration profile through time
% x = x';
% figure;
% movind = 1;
% for i = 1:size(a,1)
%     m = m0*(a(i,:).^n)./(K^n+a(i,:).^n);
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
%     plot(x,a(i,:));
%     ylim([0 10*c0]);
%     
%     drawnow
%     mov(movind,:) = getframe;
%     movind = movind+1;
% end