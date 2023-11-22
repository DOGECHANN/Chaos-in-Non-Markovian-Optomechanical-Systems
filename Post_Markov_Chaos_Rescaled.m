function dX = Post_Markov_Chaos_Rescaled(t, X, Kc, Dc, Gc, omc, P, gm)

Ac = X(1);  Bc = X(2);

f0 = Kc/2*(1-exp(-gm*t));
f1 = Kc/(2*gm)*(1-exp(-gm*t)-gm*t*exp(-gm*t));
f2 = Kc^2/(4*gm)*(1-exp(-gm*t)-gm*t*exp(-gm*t)-1/2*gm^2*t^2*exp(-gm*t));

dX = [-1i*(1+f1)*((Bc+conj(Bc))*Ac-Dc/omc*Ac+1/2)-(conj(f0)+f2)*Ac/omc; % A = X(1)
      -1i*(P/2*abs(Ac)^2+Bc)-(Gc/2)*Bc/omc % B = X(2)
      ];
end