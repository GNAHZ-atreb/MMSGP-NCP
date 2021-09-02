function [iter, cput, err, fe, suc] = m_MMSGP(tp, xk, eps, iterMAX)
%%% input
n = length(xk); sigma = 0.001; beta = 0.618; r = 0.001; alpk = ones(n,1); alpmin = 1e-2*ones(n,1);
%%% error bound
fk = fun(tp, xk); Fk = mncpfun_fb(xk, fk);
err = norm(min(xk, fk));
%%% initializations
fe = 1; iter = 0;

cput0 = tic ;
while (err > eps && iter < iterMAX)
    iter = iter + 1;
    %%% linear search
    dk = alpk.*Fk;
    zk = xk - dk;
    fzk = fun(tp, zk); Fzk = mncpfun_fb(zk, fzk); squ_dk = dk'*dk;
    squ_Fzk = Fzk'*Fzk; sqr_Fzk = sqrt(squ_Fzk); gammak = sqr_Fzk/(1+sqr_Fzk);
    left = Fzk'*dk; right = sigma*gammak*squ_dk;
    lambda = 1; m = 0 ;
    while (left < right) && (m < 50)
        m = m + 1 ;
        lambda = beta*lambda;
        zk = xk - lambda*dk; fzk = fun(tp, zk); Fzk = mncpfun_fb(zk, fzk);
        squ_Fzk = Fzk'*Fzk; sqr_Fzk = sqrt(squ_Fzk); gammak = sqr_Fzk/(1+sqr_Fzk);
        left = Fzk'*dk; right = sigma*lambda*gammak*squ_dk;
    end
    %%% new iterate
    xt = xk - ((lambda*left)/squ_Fzk)*Fzk;
    ft = fun(tp, xt); Ft = mncpfun_fb(xt, ft);
    %%% update multivarivate BB step-size vector
    sk = xt - xk; yk = Ft - Fk + r*sk;
    sk_dot_product_sk = sk.*sk; sk_dot_product_yk = sk.*yk; one_dif_sk = ones(n,1) - abs(sign(sk));
    tk = sum(sk_dot_product_sk)*one_dif_sk;  vk = sum(sk_dot_product_sk)*one_dif_sk;
    alpk = max(alpmin, (sk_dot_product_sk + tk)./(sk_dot_product_yk + vk));
    %%% update iterate
    xk = xt; Fk = Ft;
    %%% function evaluations
    fe = fe + (2+m);   
    %%% error
    err = norm(min(xt, ft));
end
%%% output
cput = toc(cput0);
if (iter == iterMAX || err > eps)
    suc = 0;
else
    suc = 1;
end