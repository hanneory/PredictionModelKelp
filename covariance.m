covAN = cov(Y_A, Y_N);
covAC = cov(Y_A, Y_C);
covNC = cov(Y_N, Y_C);
corrAN = corrcoef(Y_A, Y_N);
corrAC = corrcoef(Y_A, Y_C);
corrNC = corrcoef(Y_N, Y_C);
corrAT = corrcoef(Y_A, X_T);
corrANO3 = corrcoef(Y_A, X_XNO3);
corrAI = corrcoef(Y_A, X_I);
corrAU = corrcoef(Y_A, X_U);
corrNT = corrcoef(Y_N, X_T);
corrNNO3 = corrcoef(Y_N, X_XNO3);
corrNI = corrcoef(Y_N, X_I);
corrNU = corrcoef(Y_N, X_U);
corrCT = corrcoef(Y_C, X_T);
corrCNO3 = corrcoef(Y_C, X_XNO3);
corrCI = corrcoef(Y_C, X_I);
corrCU = corrcoef(Y_C, X_U);

covMatrix = [covAN(1,1) covAN(1,2) covAC(1,2);
             covAN(2,1) covAN(2,2) covNC(1,2);
             covAC(2,1) covNC(1,2) covNC(2,2)];

corrMatrix = ...
    [corrAN(1,1) corrAN(1,2) corrAC(1,2) corrAT(1,2) corrANO3(1,2) corrAI(1,2) corrAU(1,2);
     corrAN(2,1) corrAN(2,2) corrNC(1,2) corrNT(1,2) corrNNO3(1,2) corrNI(1,2) corrNU(1,2);
     corrAC(2,1) corrNC(2,1) corrAC(2,2) corrCT(1,2) corrCNO3(1,2) corrCI(1,2) corrCU(1,2)]

corrMatrix = ...
    [corrAN(1,1) corrAN(1,2) corrAC(1,2);
     corrAN(2,1) corrAN(2,2) corrNC(1,2);
     corrAC(2,1) corrNC(2,1) corrAC(2,2)]
