function LE = LyapunovExp(X, P, d, om, Nt)

    B1 = real(X(Nt/2:end,2));   
    format compact;
    writematrix(B1,sprintf('%g %g.txt',P, d/om));
    fname = fileread(sprintf('%g %g.txt',P, d/om));
    fclose('all');
    delete(sprintf('%g %g.txt',P, d/om));
    % datcnt = 60000;
    datcnt = length(B1);
    tau = 11;
    ndim = 4;
    ires = 10;
    maxbox = 6000;
    db = basgen(fname, tau, ndim, ires, datcnt, maxbox);
    dt = .1;
    evolve = 10;
    dismin = 0.02*range(B1);
    dismax = 0.15*range(B1);
    thmax = 30;
    [out, SUM] = fet(db, dt, evolve, dismin, dismax, thmax);
    LE = out(length(rmmissing(out))-1,4);

end