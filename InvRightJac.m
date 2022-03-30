function invJr = InvRightJac(phi)
    mag = norm(phi);
    Phi = skew(phi);
    invJr = eye(3) + 1/2 * Phi + (1/mag^2 + (1+cos(mag))/(2*mag*sin(mag))) * Phi * Phi;
end
