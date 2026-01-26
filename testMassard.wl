<<EpidCRN`;
RHSf = {-alhv*H*v + behd*H*d - muh*H, alhv*H*v - alei*i*Te - mui*i,
        lamd*d*(1-M) + lamv*v*(1-M) - mum*M, latm*M*(1-Tm) - mut*Tm - bevt*Tm*v,
        Laem*Te*M - alei*Te*i - mute*Te + bevt*Tm*v, labm*M*(1-Bm) - mub*Bm - bevb*Bm*v,
        Lapm*M*Be + lap*(1-Be) + bevb*Bm*v, -alhv*H*v + gav*i - alav*A*S*v - muv*v + C,
        -alav*A*S*v + gaa*Be - mua*A, las*(1-S)*Be, -behd*d*H + alei*i*Te + muh*H + mui*i};
varFull = {H, i, M, Tm, Te, Bm, Be, v, A, S, d};
RHS = Drop[RHSf, -1] /. d -> (1 - H - i);
var = {H, i, M, Tm, Te, Bm, Be, v, A, S};
res = EpidCRN`ODE2RN[RHS, var];
RN = res[[1]];
Print["Testing minSiph0 with Massard model..."];
result = EpidCRN`minSiph0[var, RN, RHS];
Print["SUCCESS - no infinite loop!"];
Print["mSi = ", result[[1]]];
Print["cE0 = cDFE (overdetermined system blocked)"] /; result[[2]] === result[[3]];
