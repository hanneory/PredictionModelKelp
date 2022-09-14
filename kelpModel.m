%Parameters
A_0 = 6;
alpha = 3.75*10^(-5);
C_min = 0.01;
C_struct = 0.20;
gamma = 0.5; 
epsilon = 0.22;
I_sat = 200;
J_max = 1.4*10^(-4);
k_A = 0.6;
k_dw = 0.0785;
k_C = 2.1213;
k_N = 2.72;
m_1 = 0.1085;
m_2 = 0.03;
my_max = 0.18;
N_min = 0.01;
N_max = 0.022;
N_struct = 0.01;
P_1 = 1.22*10^(-3);
P_2 = 1.44*10^(-3);
a_1 = 0.85;
a_2 = 0.3;
R_1 = 2.785*10^(-4);
R_2 = 5.429*10^(-4);
T_R1 = 285;
T_R2 = 290;
T_AP = 1694.4;
T_APH = 25924;
T_APL = 27774;
T_AR = 11033;
U_065 = 0.03;
K_X = 4;

grossPhotosynthesis(alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, 15)

function [dA, my, f_area, f_temp, f_photo, ny, dN, J, dC, P, R, E] = ...
    kelpModelFunction(A_O, alpha, C_min, C_struct, gamma, epsilon, I_sat, J_max, k_A, ...
    k_dw, k_C, k_N, m_1, m_2, my_max, N_min, N_max, N_struct, P_1, P_2, a_1, a_2, R_1, R_2, ...
    T_R1, T_R2, T_AP, T_APH, T_APL, T_AR, U_065, K_X, C, T)
    
    E = 1 - exp(gamma*(C_min-C));
    R = R_1 * exp(T_AR/T_R1-T_AR/T);


    f_temp = EffectTemp(T);
    f_area = m_1 * exp(-(A/A_O)^2)+m_2;
    my = f_area*f_photo*min([1-N_min/N, 1-C_min/C]);
    dA = (my-ny)*A;

end

function f_temp = EffectTemp (T)
    if (T>=-1.8) && (T < 10)
        f_temp = 0.08*T + 0.2;
    elseif (T>= 10) && (T<= 15)
        f_temp = 1;
    elseif (T>15) && (T<=19)
        f_temp = 19/4 - T/4;
    elseif (T>19)
        f_temp = 0;
    end
end

%function B = grossPhotosynthesis (P_1, T_AP, T_R1, T, alpha, I_sat



function P = grossPhotosynthesis ( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T)
    syms B I;
    T_PL = 271;
    T_PH = 296;
    P_max = (P_1*exp(T_AP/T_R1-T_AP/T))/(1+exp(T_AP/T- T_APL/T_PL) + exp(T_APH/T_PH - T_PH/T))

    P(B) = (alpha*I_sat) / (log(1+alpha/B)) * (alpha/(alpha+B)) * ((B/(alpha+B))^(B/alpha));
    pretty(P(B));
    Pd(B) = diff(P(B), B);
    B0 = 1*10^-9;
    nmax = 10;
    B = B0;
    for i = 1:nmax
        B = double(B - P(B)/Pd(B));
    end

    P_max2 = double(P(B))
    P_S = (alpha*I_sat)/(log(1+alpha/B));
    P = P_S * (1-exp(-(alpha*I_sat)/P_S)) * exp(- (B*I_sat)/P_S);

end



 %   while eps>=1e-5&n<=nmax

  %  for i = 1:10
   %     B(i+1) = B(i) - P(B(i))/Pd(B(i));
   %     if (abs(P(B(i+1)))<0.001)
    %        disp(double(B(B+1)));
     %       break;
      %  end
   % end



