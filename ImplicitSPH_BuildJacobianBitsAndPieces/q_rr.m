function out = q_rr(dist3, h)
q = norm(dist3,2) / h;
out = zeros(3);
out = 1.0 / q * (1.0 / h^2 * eye(3) - 1.0 / (h^4 * q^2) * dist3 * dist3');
