function dX = Post_Markov(t, X, d, P, Kappa, Gamma, om, gm1) %Post_Markov(t, X, d, P, Kappa, Gamma, om, g0, gm1)

Ac = X(1);  Bc = X(2);

f0 = Kappa/2*(1-exp(-gm1*t));
f1 = Kappa/(2*gm1)*(1-exp(-gm1*t)-gm1*t*exp(-gm1*t));
f2 = Kappa^2/(4*gm1)*(1-exp(-gm1*t)-gm1*t*exp(-gm1*t)-1/2*gm1^2*t^2*exp(-gm1*t));

dX = [-1i*(1+f1)*((Bc+conj(Bc))*Ac-d/om*Ac+1/2)-(conj(f0)+f2)*Ac/om; % A = X(1)
      -1i*(P/2*abs(Ac)^2+Bc)-(Gamma/2)*Bc/om % B = X(2)
      ];
end