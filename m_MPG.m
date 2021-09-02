function [iter, cput, err, fe, suc] = m_MPG(tp, xk, eps, iterMAX)
%%% input
delta = 1.01;
lambdak = 0.4; alpha = 0.4/delta;
%%% error bound
fk = fun(tp, xk);
err = norm(min(xk, fk));

%%% initializations
yk = xk; fyk = fk;
iter = 0; fe = 1;

%%% the start of time
cput0 = tic;
while (err > eps && iter < iterMAX)
    iter = iter + 1;
    %%% new ietrate
    dk = - min(xk, lambdak*fyk);
    xt = xk + dk;
    
    yt = xt + delta*dk; fyt = fun(tp, yt);
    %%% parameter update
    sqr_df = norm(fyt - fyk);
    if (sqr_df > 1e-6)
        lambdak = min(lambdak, alpha*norm(yt - yk)/sqr_df);
    end
    %%% update
    xk = xt; fk = fun(tp, xt);
    fyk = fyt;
    
    %%% err bound
    err = norm(min(xk, fk));
    %%% the number of function evaluations
    fe = fe + 2;
end
%%% output
cput = toc(cput0);
if (iter == iterMAX || err > eps)
    suc = 0;
else
    suc = 1;
end