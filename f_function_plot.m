om=10^0;%Unit of measuring time
Kappa=1*om;%Optical damping
Gamma=10^(-3)*om;%Mechanical damping
g0=0.1*om; % coupling constant

gm1 = 2*om; % non-Markovian correlation frequency \gamma

Nt = 100; % Number of points for time
t_0 = 0; % Initial time
tspan = 5; % Final time
t = linspace(t_0, tspan, Nt); % Time vector

f0 = Kappa/2*(1-exp(-gm1.*t));
f1 = Kappa/(2*gm1)*(1-exp(-gm1.*t)-gm1.*t.*exp(-gm1.*t));
f2 = Kappa^2/(4*gm1)*(1-exp(-gm1.*t)-gm1.*t.*exp(-gm1.*t)-1/2*gm1^2.*t.^2.*exp(-gm1.*t));

figure 

plot(t(1:Nt),f0,'LineWidth',1);
hold on
plot(t(1:Nt),f1,'LineWidth',1);
hold on
plot(t(1:Nt),f2,'LineWidth',1);
hold off

legend('f0','f1','f2');
axis([t_0 tspan 0 0.6]);
set(findall(gcf,'-property','FontSize'),'FontSize',14);