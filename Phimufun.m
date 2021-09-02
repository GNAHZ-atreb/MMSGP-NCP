function phi = Phimufun(x, f, mu)

st = (x-f).^2 + 4*(mu^2);
phi = x + f - sqrt(st);

end