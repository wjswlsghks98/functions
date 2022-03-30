function R = Exp_map(phi)
%EXP Returns exponential map from R3 --> SO3 
    Phi = skew(phi);
    mag = norm(phi);
    if mag == 0
        R = eye(3);
    else
        R = eye(3) + sin(mag)/mag * Phi + (1-cos(mag))/mag^2 * Phi * Phi;
    end
end
