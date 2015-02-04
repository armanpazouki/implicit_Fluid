function out = dW(dist3,h)
q = norm(dist3,2) / h;
if (q < 1)
    out = -3*(2-q)^2 + 12*(1-q)^2;
elseif (q < 2)
    out = -3*(2-q)^2;
else
    out = 0;
end
out = out * (1.0 / (4 * pi * h^3));