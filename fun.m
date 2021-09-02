function f = fun(tp, x)
n =  length(x);
%%%
switch tp
    case{1}
        f = exp(x) - 1;
    case{2}
        f = atan(x) - pi/6;
    case{3}
        ps = n\(1:n)';
        f = exp(ps.*x) - 1;
    case{4}
        f = log(x + 1) - n\x;
    case{5}
        f = zeros(4,1);
        sx1 = x(1)^2; sx2 = x(2)^2; 
        f(1) = 3*sx1 + 2*x(1)*x(2) + 2*sx2 + x(3) + 3*x(4) - 6;
        f(2) = 2*sx1 + x(1) + sx2 + 3*x(3) + 2*x(4) - 2;
        f(3) = 3*sx1 + x(1)*x(2) + 2*sx2 + 2*x(3) + 3*x(4) - 1;
        f(4) = sx1 + 3*sx2 + 2*x(3) + 3*x(4) - 3;
    case{6}
        f = 2*x - sin(abs(x));
    case{7}
        f = x  - sin(abs(x-1));
    case{8}
        f = log(abs(x) + 1) - n\x;
    case{9}
        xt1 = [0; x(1:n-1)]; xt2 = [x(2:n); 0]; h = 1/(n+1);
        f = x - exp(cos(h*(xt1 + x + xt2)));
    case{10}
        xt1 = [0; x(1:n-1)]; xt2 = [x(2:n); 0]; h = 1/(n+1);
        f = x - exp(cos(h*abs(xt1 + x + xt2)));
    case{11}
        abs_x = abs(x);
        f = min(min(abs_x, x.^2), max(abs_x, x.^3));
    case{12}
        x_squ = x.^2; xt1 = [0; x(1:n-1)]; xt2 = [x_squ(2:n); 0];
        t = (n+1)\(1:n)';
        f = t.*x.^3 - n\abs(xt1) + n\abs(xt2);
    case{13}
        f = exp(abs(x)) - ones(n,1);
    case{14}
        ct = 1./(1:n)'; xt1 = [x(2:n); 0]; xt2 = [0; x(1:n-1)];
        f = x + atan(x) - cos(ct.*(xt1 + xt2)) - ones(n,1);
    case{15}
        ct = n\(1:n)'; xt1 = exp(ct.*x); sinhx = 0.5*(xt1 - 1./xt1);
        xt2 = [0; x(1:n-1)]; xt3 = [x(2:n); 0];
        f = sinhx + xt2 + xt3 + ct;
    otherwise
        disp( 'This input index "tp" does not exist in this test library. Please re-enter£¡' );
end
