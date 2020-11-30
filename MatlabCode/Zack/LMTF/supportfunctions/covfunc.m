function f = covfunc()
f = {@covProd,{{@covMask,{[1 0],{@covMaterniso,1}}},{@covMask,{[0 1],@covPeriodPi}}}};
