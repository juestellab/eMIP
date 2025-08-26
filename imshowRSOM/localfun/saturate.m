function th = saturate(th, ContrastFactorLow, ContrastFactorHigh, ContrastFactorMethod)
    switch ContrastFactorMethod
    case 'relative'
        f = th(2)-th(1);
        th(2) = ContrastFactorHigh*(f)+th(1);
        th(1) = ContrastFactorLow*(f)+th(1);
    case 'absolute'
        th(1) = ContrastFactorLow*th(1);
        th(2) = ContrastFactorHigh*th(2);
    end
end
