tic

clear all

om=10^0;%Unit of measuring time
Kappa=1*om;%Optical damping
Gamma=10^(-3)*om;%Mechanical damping
g0=0.1*om; % coupling constant

Nt = 20000; % Number of points for time
t_0 = 0; % Initial time
tspan = 2000; % Final time
t = linspace(t_0, tspan/om, Nt); % Time vector

% p = 1.5;
% aL=((om^2)/(2*g0))*sqrt(p/2);
% aL = 4.15;
P = 1.4;

gridsize = 1000;
Dleft = -1.14;
Dright = -0.4;
LinD = linspace(Dleft, Dright, gridsize);
% gammaline = [10 2];
gm1 = 2; % non-Markovian correlation frequency \gamma


Lm = zeros(gridsize);
Dm = zeros(gridsize);

Lmm = zeros(gridsize);
Dmm = zeros(gridsize);

figure
tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'none');

% for jj = 1:2
%     gm1 = gammaline(jj)*om; %Optical memory frequency
    
nexttile 

for ii = 1:gridsize
    d = LinD(ii);
    [t,X]=ode45(@(t,X) Post_Markov(t, X, d, P, Kappa, Gamma, om, gm1), t, [0 0]);
%%
    LCE = LyapunovExp(X, P, d, om, Nt);
    Lm(ii) = LCE;
    Dm(ii) = d;

%%
    X2 = abs(X(Nt/2:end,1));   
    peakvalues = findpeaks(X2);
    valleyvalues = -findpeaks(-X2);
    stableB = [peakvalues; valleyvalues]; 
    lp = length(stableB);
    y = stableB;
    z = ones(lp,1);
    plot(d*z, y, '.', 'markersize', 1, 'color', 'k');
    hold on
end
hold off
set(gca,'XTick',[]);
% title(sprintf('\\gamma = %g', gm1/om));
%xlabel('\Delta');
ylabel('Amplitude');

axis ([Dleft Dright 0 1])

nexttile
plot(Dm, Lm, '-');
%title(sprintf('LE Post Markov P = %g \\gamma = %g', [p gm1]));
xlabel('\Delta');
ylabel('LE');

axis ([Dleft Dright -0.05 0.2])
set(findall(gcf,'-property','FontSize'),'FontSize',14);

% end
toc
