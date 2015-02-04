function out = Distance(r_a, r_b)
global cMin cMax
side3 = cMax - cMin;
out = r_a - r_b;
if out(1) > 0.5 * side3(1)
    out(1) = out(1) - side3(1);
end
if out(2) > 0.5 * side3(2)
    out(2) = out(2) - side3(2);
end
if out(3) > 0.5 * side3(3)
    out(3) = out(3) - side3(3);
end

if out(1) < -0.5 * side3(1)
    out(1) = out(1) + side3(1);
end
if out(2) < -0.5 * side3(2)
    out(2) = out(2) + side3(2);
end
if out(3) < -0.5 * side3(3)
    out(3) = out(3) + side3(3);
end
