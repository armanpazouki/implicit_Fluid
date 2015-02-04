function out = Perturb(vec3, thresh)
out = vec3 + thresh.*[2 * rand() - 1; 2 * rand() - 1; 2 * rand() - 1];
