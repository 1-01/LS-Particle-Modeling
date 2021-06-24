 function dt = calc_dt(Ubar, sigu, tau)
    delx = 1;
    DT_FACT1 = 0.025;
    Courant = DT_FACT1*tau*(Ubar+sigu)/delx;
    if Courant < 0.5
        DT_FACT = DT_FACT1;
    else
        DT_FACT = 0.5*delx/(tau*(Ubar + sigu));
    end
    dt = DT_FACT*tau;
end