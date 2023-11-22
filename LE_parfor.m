tic

%%
% Setup

omc = 10^0; % Unit of measuring time
Kc = 1*10^(0)*omc; % Optical damping
Gc = 1*10^(-3)*omc; % Mechanical damping

Nt = 20000; % Number of points for time
t_0 = 0; % Initial time
tspan = 2000; % Final time

gridsize = 200;
LinD = linspace(-1.4, -0.4, gridsize);
LinP = linspace(0.8, 1.6, gridsize);
% LinAL = ((om^2)/(2*g0))*sqrt(LinP./2);
[M, N] = meshgrid(LinD, LinP);

gm = 10; % non-Markovian correlation frequency \gamma

figure

nexttile;
Z1 = nan(size(N));
c1 = surf(M, N, Z1,'EdgeColor','none');
colorbar
clim([0 0.14]);
colormap(turbo);
ylabel('P');
xlabel('\Delta'); 
view(0, 90)
axis([-1.4 -0.4 0.8 1.6]);
% title(sprintf('LE, \\gamma = %g', gm1/om));
D1 = parallel.pool.DataQueue;
D1.afterEach(@(y) updateSurface(c1, y));
set(findall(gcf,'-property','FontSize'),'FontSize',14);

% hold on % gamma = 10
% plot(-0.89,0.8,'y*'); 
% plot(-0.95,0.9,'y*'); 
% plot(-1.00,1.0,'y*'); 
% plot(-1.05,1.1,'y*');  
% plot(-1.10,1.2,'y*');  
% plot(-1.14,1.3,'y*'); 

% plot(-1.18,1.4,'w*');  
% plot(-1.22,1.5,'w*');  
% plot(-1.26,1.6,'w*');
% 
% plot(-0.69,1.08,'w*'); plot(-0.61,1.08,'w*');
% plot(-0.77,1.1,'w*'); plot(-0.57,1.1,'w*');
% plot(-0.93,1.2,'w*'); plot(-0.50,1.2,'w*');
% plot(-1.08,1.3,'w*'); plot(-0.48,1.3,'w*');
% plot(-0.46,1.4,'w*');
% plot(-0.46,1.5,'w*');
% plot(-0.45,1.6,'w*');
% hold off


% hold on % gamma = 2
% % plot(-0.89,0.8,'y*'); 
% % plot(-0.95,0.9,'y*'); 
% % plot(-1.00,1.0,'y*');  
% % plot(-1.05,1.1,'y*');
% % plot(-1.10,1.2,'y*');
% 
% plot(-1.14,1.3,'w*');  
% plot(-1.18,1.4,'w*');
% plot(-1.22,1.5,'w*');
% plot(-1.25,1.6,'w*');
% 
% plot(-0.77,1.0,'w*'); plot(-0.48,1.0,'w*');  
% plot(-0.91,1.1,'w*'); plot(-0.46,1.1,'w*');  
% plot(-1.04,1.2,'w*'); plot(-0.44,1.2,'w*');
% plot(-0.44,1.3,'w*');  
% plot(-0.45,1.4,'w*');  
% plot(-0.46,1.5,'w*');  
% plot(-0.47,1.6,'w*');
% plot(-0.69,0.96,'w*'); plot(-0.51,0.96,'w*');
% plot(-0.61,0.935,'w*'); plot(-0.55,0.935,'w*');
% hold off

parfor ii = 1:numel(N)   
    Dc = M(ii);
    P = N(ii);
%     aL=((om^2)/(2*g0))*sqrt(p/2);
%     P = 8*aL^2*g0^2/om^4
    t = linspace(t_0, tspan, Nt); % Time vector   
    [t,X]=ode45(@(t,X) Post_Markov_Chaos_Rescaled(t, X, Kc, Dc, Gc, omc, P, gm), t, [0 0]);
    LE = LyapunovExp (X, P, Dc, omc, Nt);
    send(D1, [ii LE]);
end

toc