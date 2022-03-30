function phi = Log_map(R)
%EXP Returns exponential map from SO3 --> R3 
    psi = acos((trace(R)-1)/2);
    logR = psi * (R - R')/(2*sin(psi));
    phi = [logR(3,2); logR(1,3); logR(2,1)];
end
