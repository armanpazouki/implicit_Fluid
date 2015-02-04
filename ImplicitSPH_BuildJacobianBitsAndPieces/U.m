function out = U(r, v, rpm_a, rpm_b)
global etha
global h;
global m;
rho_a = rpm_a(1);
p_a = rpm_a(2);
mu_a = rpm_a(3);
rho_b = rpm_b(1);
p_b = rpm_b(2);
mu_b = rpm_b(3);
 
denomSquare = (norm(r))^2 + etha^2;
dw = dW(r, h);
qr = q_r(r, h);
gradW = dw * qr;
out = m * ( p_a / (rho_a^2) + p_b / (rho_b^2) ) * gradW' ...
    - m * (mu_a + mu_b) / (0.5*(rho_a + rho_b))^2 * 1.0 / denomSquare * (r' * gradW') * v;
