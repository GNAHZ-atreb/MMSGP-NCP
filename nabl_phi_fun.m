function [nap, nbp] = nabl_phi_fun(x, f, lambda, p)
n = length(x); nap = zeros(n, 1); nbp = zeros(n, 1);

ax = abs(x); af = abs(f);
xfp = ax.^p + af.^p; xfpt = xfp.^(1 - 1/p);

for i = 1:n
    if (xfpt(i) > 1e-6)
        phi_pi = x(i) + f(i) - xfp(i)^(1/p);
        x_plusi = max(x(i), 0); f_plusi = max(f(i), 0);
        phi_lamb_pi = lambda*phi_pi + (1 - lambda)*(x_plusi*f_plusi);
        
        axti = ax(i)^(p-1);
        ta1i = 1 - (sign(x(i)) * axti / xfpt(i)); ta2i = f_plusi.*sign(ax(i));
        nap(i) = phi_lamb_pi * (lambda*ta1i + (1-lambda)*ta2i);

        afti = af(i)^(p-1);
        tb1i = 1 - (sign(f(i)) * afti / xfpt(i)); tb2i = x_plusi.*sign(af(i));
        nbp(i) = phi_lamb_pi * (lambda*tb1i + (1-lambda)*tb2i);
    end
end

