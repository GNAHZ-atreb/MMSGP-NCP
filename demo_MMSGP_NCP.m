clc;clear;warning('off')
%%%%%%%%%%%%%%%%%% initializations
tp = 12; n = 20000; ip_num = 100;
eps = 1e-4; iterMAX = 10000;
%%% initial points
if (tp == 5)
    n = 4;
end
%%% initial points
X0 = rand(n,ip_num);
%%% 
iterMSGP = zeros(ip_num,1); cputMSGP = iterMSGP; errMSGP = iterMSGP; feMSGP = iterMSGP; sucMSGP = iterMSGP;
iterNDF = iterMSGP; cputNDF = iterMSGP; errNDF = iterMSGP; feNDF = iterMSGP; sucNDF = iterMSGP;
iterMPG = iterMSGP; cputMPG = iterMSGP; errMPG = iterMSGP; feMPG = iterMSGP; sucMPG = iterMSGP;
iterMMSGP = iterMSGP; cputMMSGP = iterMSGP; errMMSGP = iterMSGP; feMMSGP = iterMSGP; sucMMSGP = iterMSGP;

for ip = 1:ip_num
    x0 = X0(:,ip);
    %save x0_1; save x0_2; save x0_3; save x0_4; save x0_5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% methods & output
    fprintf('method        iter       cput        err           fe         s/f \n')
    %%% MSGP
%     [iterMSGP(ip), cputMSGP(ip), errMSGP(ip), feMSGP(ip), sucMSGP(ip)] = m_MSGP(tp, x0, eps, iterMAX);
%     fprintf('MSGP          %-6d     %-7.4f     %3.1e       %-8d    %-6d \n', iterMSGP(ip), cputMSGP(ip), errMSGP(ip), feMSGP(ip), sucMSGP(ip))
    %%% NDF
%     [iterNDF(ip), cputNDF(ip), errNDF(ip), feNDF(ip), sucNDF(ip)] = m_NDF(tp, x0, eps, iterMAX);
%     fprintf('NDF           %-6d     %-7.4f     %3.1e       %-8d    %-6d \n', iterNDF(ip), cputNDF(ip), errNDF(ip), feNDF(ip), sucNDF(ip))
    %%% MPG
%     [iterMPG(ip), cputMPG(ip), errMPG(ip), feMPG(ip), sucMPG(ip)] = m_MPG(tp, x0, eps, iterMAX);
%     fprintf('MPG           %-6d     %-7.4f     %3.1e       %-8d    %-6d \n', iterMPG(ip), cputMPG(ip), errMPG(ip), feMPG(ip), sucMPG(ip))
    %%% MMSGP
    [iterMMSGP(ip), cputMMSGP(ip), errMMSGP(ip), feMMSGP(ip), sucMMSGP(ip)] = m_MMSGP(tp, x0, eps, iterMAX);
    fprintf('MMSGP         %-6d     %-7.4f     %3.1e       %-8d    %-6d \n', iterMMSGP(ip), cputMMSGP(ip), errMMSGP(ip), feMMSGP(ip), sucMMSGP(ip))

    fprintf('------------------------------------------------------------------------ ip%-2d \n', ip )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% methods & avg_output
fprintf('method         avg_iter      avg_cput       avg_err         avg_fe        suc_rate   \n')
%%% MSGP
[avg_iterMSGP, avg_cputMSGP, avg_errMSGP, avg_feMSGP]  = valid_mean_fun(iterMSGP, cputMSGP, errMSGP, feMSGP, sucMSGP);
fprintf('MSGP           %-6.1f        %-7.4f        %3.1e         %-8.1f      %.2f%% \n', avg_iterMSGP, avg_cputMSGP, avg_errMSGP, avg_feMSGP, sum(sucMSGP)/ip_num*100)
%%% NDF
[avg_iterNDF, avg_cputNDF, avg_errNDF, avg_feNDF]  = valid_mean_fun(iterNDF, cputNDF, errNDF, feNDF, sucNDF);
fprintf('NDF            %-6.1f        %-7.4f        %3.1e         %-8.1f      %.2f%% \n', avg_iterNDF, avg_cputNDF, avg_errNDF, avg_feNDF, sum(sucNDF)/ip_num*100)
%%% MPG
[avg_iterMPG, avg_cputMPG, avg_errMPG, avg_feMPG]  = valid_mean_fun(iterMPG, cputMPG, errMPG, feMPG, sucMPG);
fprintf('MPG            %-6.1f        %-7.4f        %3.1e         %-8.1f      %.2f%% \n', avg_iterMPG, avg_cputMPG, avg_errMPG, avg_feMPG, sum(sucMPG)/ip_num*100)
%%% MMSGP
[avg_iterMMSGP, avg_cputMMSGP, avg_errMMSGP, avg_feMMSGP]  = valid_mean_fun(iterMMSGP, cputMMSGP, errMMSGP, feMMSGP, sucMMSGP);
fprintf('MMSGP          %-6.1f        %-7.4f        %3.1e         %-8.1f      %.2f%% \n', avg_iterMMSGP, avg_cputMMSGP, avg_errMMSGP, avg_feMMSGP, sum(sucMMSGP)/ip_num*100)
fprintf('--------------------------------------------------------------------------------- avg_output \n' )