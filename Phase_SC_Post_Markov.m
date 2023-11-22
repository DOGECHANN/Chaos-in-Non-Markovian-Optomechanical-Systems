clear all

om=10^0;%Unit of measuring time
Kappa=1*om;%Optical damping
Gamma=10^(-3)*om;%Mechanical damping
g0=0.1*om; % coupling constant


d = -0.65;

P = 1.37;

KB = 1.38e-23;
T = 0.002;
hbar = 1.05457e-34;

Nt = 100000; % Number of points for time
t_0 = 0; % Initial time
tspan = 10000; % Final time

% gammaline = [10 2 1 0.8];
% 
figure
tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'compact');
% for jj = 1:4
% 
%     gm1 = gammaline(jj)*om;

gm1 = 2*om; % non-Markovian correlation frequency \gamma


%     aL=((om^2)/(2*g0))*sqrt(p/2);
    t = linspace(t_0, tspan, Nt); % Time vector
    
    [t,X]=ode45(@(t,X) Post_Markov(t, X, d, P, Kappa, Gamma, om, gm1), t, [0 0]);

   
    A1 = real(X(Nt/2:end,1));
    A2 = imag(X(Nt/2:end,1));
    B1 = real(X(Nt/2:end,2));
    B2 = imag(X(Nt/2:end,2));
    
    
nexttile
    patch(A1, A2, B1, B2,'edgecolor','flat','facecolor','none','Marker','.','markersize',5,'LineStyle','none');
    xlabel('Z_1')
    ylabel('Z_2')
    zlabel('Z_3')
    cb = colorbar;                                  % create and label the colorbar
%     cb.Label.String = 'Im \langle b\rangle';
%     title(sprintf('\\gamma = %g', gm1/om));
    view(30, 30)
   
% end

LCE = LyapunovExp(X, P, d, om, Nt);
set(findall(gcf,'-property','FontSize'),'FontSize',14);

Nt = 100; % Number of points for time
t_0 = 0; % Initial time
tspan = 5; % Final time
t = linspace(t_0, tspan, Nt); % Time vector

f0 = Kappa/2*(1-exp(-gm1.*t));
f1 = Kappa/(2*gm1)*(1-exp(-gm1.*t)-gm1.*t.*exp(-gm1.*t));
f2 = Kappa^2/(4*gm1)*(1-exp(-gm1.*t)-gm1.*t.*exp(-gm1.*t)-1/2*gm1^2.*t.^2.*exp(-gm1.*t));

nexttile

plot(t(1:Nt),f0,'LineWidth',1);
hold on
plot(t(1:Nt),f1,'LineWidth',1);
hold on
plot(t(1:Nt),f2,'LineWidth',1);
hold off

legend('f_0','f_1','f_2');
axis([t_0 tspan 0 0.6]);
set(findall(gcf,'-property','FontSize'),'FontSize',14);