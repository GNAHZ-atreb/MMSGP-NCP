function [iter, cput, err, fe, suc] = m_NDF(tp, xk, eps, iterMAX)
%%% input
n = length(xk);
sigma = 1.0e-10; ms = 5; M = 5; gamma = 0.1;
beta = 0.2; lambda = 0.9; p = 2;
%%% error bound
fk = fun(tp, xk); Philpk = Phi_lamb_p_fun(xk, fk, lambda, p); squ_Philpk = 0.5*(Philpk'*Philpk);
err = norm(min(xk, fk));
%%% initializations
fe = 1; iter = 0; ctm = zeros(n ,1);

cput0 = tic ;
while (err > eps && iter < iterMAX)
    iter = iter + 1; ctm(iter) = squ_Philpk;
    %%% maxmum of merit function from k-th iteration to (k-mk)-th iteration  
    if (iter <= ms)
        mk = 0;
    else
        mk = min(M, mk+1); 
    end
    squ_ell_Philpk = max(ctm(iter-mk: iter));
    %%% search direction
    [da, db] = nabl_phi_fun(xk, fk, lambda, p);
    sum_dab = da + db; squ_sum_dab = sum_dab'*sum_dab;
    dk = -sum_dab;
    %%% initial update
    xt = xk + dk; ft = fun(tp, xt);
    Philpt = Phi_lamb_p_fun(xt, ft, lambda, p); squ_Philpt = 0.5*(Philpt'*Philpt);
    %%% linear search
    right = squ_ell_Philpk - sigma*squ_sum_dab;
    intk = 0;
    while (squ_Philpt > right && intk < 50)
        intk = intk + 1; gamma_s = gamma^intk; beta_s = beta^intk;
        dk = -db - beta_s*da;
        xt = xk + gamma_s*dk; ft = fun(tp, xt);
        Philpt = Phi_lamb_p_fun(xt, ft, lambda, p); squ_Philpt = 0.5*(Philpt'*Philpt);
    
        right = squ_ell_Philpk - sigma*(gamma_s^2)*squ_sum_dab;
    end
    %%% update iterate
    xk = xt; fk = ft; squ_Philpk = squ_Philpt;
    %%% function evaluations
    fe = fe + (1+intk);   
    %%% error
    err = norm(min(xk, fk));
end
%%% output
cput = toc(cput0);
if (iter == iterMAX || err > eps)
    suc = 0;
else
    suc = 1;
end
