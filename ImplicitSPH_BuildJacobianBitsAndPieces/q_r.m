function out = q_r(dist3, h)
q = norm(dist3,2) / h;
out = zeros(1,3);
out = (1.0 / (h^2 * q)) * dist3';