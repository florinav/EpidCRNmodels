<<EpidCRN`;
RHS = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
var = {s, i};
res = EpidCRN`ODE2RN[RHS, var];
RN = res[[1]];
rts = res[[2]];
Print["Testing minSiph0:"];
result = EpidCRN`minSiph0[var, RN, RHS];
Print["result=", result];
