function [iter, cput, err, fe, suc] = m_MSGP(tp, xk, eps, iterMAX)
%%% input
n = length(xk); sigma = 0.001; epsilon = 1e-10; invepsilon = 1e+10; 
rho = 0.5; r = 0.01; tau = 1e-8; alpk = ones(n,1);
%%% error bound
fk = fun(tp, xk); Fk = mncpfun_fb(xk, fk);
err = norm(min(xk, fk));
%%% initializations
n = length(xk);
fe = 1; iter = 0;

cput0 = tic ;
while (err > eps && iter < iterMAX)
    iter = iter + 1;
    %%% linear search
    dk = alpk.*Fk; squ_dk = dk'*dk;
    zk = xk - dk;
    fzk = fun(tp, zk); Fzk = mncpfun_fb(zk, fzk); squ_Fzk = Fzk'*Fzk; sqr_Fzk = sqrt(squ_Fzk);
    left = Fzk'*dk; right = sigma*sqr_Fzk*squ_dk;
    %%% compute the trial steplength
    tauk = xk - tau*dk; ftauk = fun(tp, tauk); Ftauk = mncpfun_fb(tauk, ftauk);
    beta = ((Fk'*dk)/(dk'*(Ftauk - Fk)))/tau;
    if (beta <= 1e-6)
        beta = 1;
    end
    lambda = 1; m = 0 ;
    while (left < right && m < 50)
        m = m + 1 ;
        lambda = rho*lambda;
        zk = xk - (beta*lambda)*dk;
        fzk = fun(tp, zk); Fzk = mncpfun_fb(zk, fzk); squ_Fzk = Fzk'*Fzk; sqr_Fzk = sqrt(squ_Fzk);
        left = Fzk'*dk; right = sigma*beta*lambda*sqr_Fzk*squ_dk;
    end
    %%% new iterate
    xt = xk - ((lambda*left)/squ_Fzk)*Fzk;
    ft = fun(tp, xt); Ft = mncpfun_fb(xt, ft);
    %%% update BB vector
    sk = xt - xk; yk = Ft - Fk + r*sk; bt = (sk'*sk)/(sk'*yk);
    for i = 1:n
        if (yk(i)*sk(i) > 0)
            alpk(i) = yk(i)/sk(i);
        else
            alpk(i) = bt;
        end
        if (alpk(i) <= epsilon || alpk(i) >= invepsilon)%% delta
            if (err > 1)
                alpk(i) = 1;
            elseif (err < 1e-5)
                alpk(i) = 1e+5;
            else
                alpk(i) = 1/err;
            end
        end
    end
    %%% uodate iterate
    xk = xt; Fk = Ft;
    %%% function evaluations
    fe = fe + (3+m);   
    %%% error
    err = norm(min(xt, ft));
end
%%% output
cput = toc(cput0);
if (iter == iterMAX || err > eps || ismissing(err) == 1)
    suc = 0;        
else
    suc = 1;
end