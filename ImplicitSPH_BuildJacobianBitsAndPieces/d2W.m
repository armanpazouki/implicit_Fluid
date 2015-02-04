function out = d2W(dist3,h)
q = norm(dist3,2) / h;
if (q < 1)
    out = 6*(2-q) - 24*(1-q);
elseif (q < 2)
    out = 6*(2-q);
else
    out = 0;
end
out = out * (1.0 / (4 * pi * h^3));