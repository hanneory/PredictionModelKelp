A = Y_A(:, 9000);
N = Y_N(:, 9000);
C = Y_C(:, 9000);
T = X_T(:, 9000);
NO3 = X_XNO3(:, 9000);
I = X_I(:,9000);
U = X_U(:,9000);

covAN = cov(A, N);
covAC = cov(A, C);
covNC = cov(N, C); 
covAT = cov(A, T);
covANO3 = cov(A, NO3);
covAI = cov(A, I);
covAU = cov(A, U);
covNT = cov(N,T);
covNNO3 = cov(N,NO3);
covNI = cov(N,I);
covNU = cov(N,U);
covCT = cov(C,T);
covCNO3 = cov(C,NO3);
covCI = cov(C,I);
covCU = cov(C,U);
covTNO3 = cov(T, NO3);
covTI = cov(T, I);
covTU = cov(T, U);
covNO3I = cov(NO3, I);
covNO3U = cov(NO3, U);
covIU = cov(I,U);

corrAN = corrcoef(A, N);
corrAC = corrcoef(A, C);
corrNC = corrcoef(N, C);
corrAT = corrcoef(A, T);
corrANO3 = corrcoef(A, NO3);
corrAI = corrcoef(A, I);
corrAU = corrcoef(A, U);
corrNT = corrcoef(N, T);
corrNNO3 = corrcoef(N, NO3);
corrNI = corrcoef(N, I);
corrNU = corrcoef(N, U);
corrCT = corrcoef(C, T);
corrCNO3 = corrcoef(C, NO3);
corrCI = corrcoef(C, I);
corrCU = corrcoef(C, U);
corrTNO3 = corrcoef(T, NO3);
corrTI = corrcoef(T, I);
corrTU = corrcoef(T, U);
corrNO3I = corrcoef(NO3, I);
corrNO3U = corrcoef(NO3, U);
corrIU = corrcoef(I, U);

covMatrix = [covAN(1,1) covAN(1,2) covAC(1,2) covAT(1,2) covANO3(1,2) covAI(1,2) covAU(1,2);
             covAN(2,1) covAN(2,2) covNC(1,2) covNT(1,2) covNNO3(1,2) covNI(1,2) covNU(1,2);
             covAC(2,1) covNC(1,2) covNC(2,2) covCT(1,2) covCNO3(1,2) covCI(1,2) covCU(1,2);
             covAT(2,1) covNT(2,1) covCT(2,1) covAT(2,2) covTNO3(1,2) covTI(1,2) covTU(1,2);
             covANO3(2,1) covNNO3(2,1) covCNO3(2,1) covTNO3(2,1) covANO3(2,2) covNO3I(1,2) covNO3U(1,2);
             covAI(2,1) covNI(2,1) covCI(2,1) covTI(2,1) covNO3I(2,1) covAI(2,2) covIU(1,2);
             covAU(2,1) covNU(2,1) covCU(2,1) covTU(2,1) covNO3U(2,1) covIU(2,1) covAU(2,2)];

corrMatrix = ...
    [corrAN(1,1) corrAN(1,2) corrAC(1,2) corrAT(1,2) corrANO3(1,2) corrAI(1,2) corrAU(1,2);
     corrAN(2,1) corrAN(2,2) corrNC(1,2) corrNT(1,2) corrNNO3(1,2) corrNI(1,2) corrNU(1,2);
     corrAC(2,1) corrNC(2,1) corrAC(2,2) corrCT(1,2) corrCNO3(1,2) corrCI(1,2) corrCU(1,2);
     corrAT(2,1) corrNT(2,1) corrCT(2,1) corrAT(2,2) corrTNO3(1,2) corrTI(1,2) corrTU(1,2);
     corrANO3(2,1) corrNNO3(2,1) corrCNO3(2,1) corrTNO3(2,1) corrANO3(2,2) corrNO3I(1,2) corrNO3U(1,2)
     corrAI(2,1) corrNI(2,1) corrCI(2,1) corrTI(2,1) corrNO3I(2,1) corrAI(2,2) corrIU(1,2);
     corrAU(2,1) corrNU(2,1) corrCU(2,1) corrTU(2,1) corrNO3U(2,1) corrIU(2,1) corrAU(2,2)];

corrMatrixShort = ...
    [corrAN(1,1) corrAN(1,2) corrAC(1,2);
     corrAN(2,1) corrAN(2,2) corrNC(1,2);
     corrAC(2,1) corrNC(2,1) corrAC(2,2)];
