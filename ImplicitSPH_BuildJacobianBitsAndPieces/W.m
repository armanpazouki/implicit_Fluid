function out = W(dist3,h)
q = norm(dist3,2) / h;
out = 0;
if (q < 1)
    out = (2-q)^3 - 4*(1-q)^3;
elseif (q < 2)
    out = (2-q)^3;
else
    out = 0;
end
out = out * (1.0 / (4 * pi * h^3));
