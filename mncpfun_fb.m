function F = mncpfun_fb(x, f)

F = f + x - sqrt(x.^2 + f.^2);

end