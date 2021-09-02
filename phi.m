function phit = phi( tp, x )

fx = fun( tp,x );
fbnft = sqrt(x.^2 + fx.^2) - (x+fx);
phit = 0.5* ( fbnft'*fbnft );

end