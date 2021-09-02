function phi_NDF = Phi_lamb_p_fun(x, f, lambda, p)

ax = abs(x); af = abs(f);
phi_p = x + f - (ax.^p + af.^p).^(1/p);
x_plus = max(x, 0); f_plus = max(f, 0);

phi_NDF = lambda*phi_p + (1-lambda)*(x_plus.*f_plus);

end