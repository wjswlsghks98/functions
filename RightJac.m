function Jr = RightJac(phi)
%JR is the Right Jacobian of SO(3)  
    mag = norm(phi);
    Phi = skew(phi);
    Jr = eye(3) - (1-cos(mag))/mag^2 * Phi + (mag - sin(mag))/mag^3 * Phi * Phi;
end
