% Analytical Poiseuille Profile based on series solution
function out = V(z, t)
global F L nu
% F: body force
% L: Channel width
% nu: kinematic viscosity
% z: distance from lower wall (vector)
% t: time

out = -F / (2*nu)  * z .* (z - L);
for n=0:10
    out = out - (4 * F * L ^ 2 / (nu * pi^3 * (2 * n + 1)^3)) * ...
        sin(pi * z * (2 * n + 1) / L) * ...
        exp(-(2 * n + 1) ^ 2 * pi ^ 2 * nu / L ^ 2 * t);
end
