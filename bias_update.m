function bias_n = bias_update(bias, dt, str)
    %% Bias Update Function: Brownian Motion
    % reference: IMU Preintegration on Manifold for Efficient
    % Visual-Inertial Maximum-a-Posteriori Estimation
    
    % Std acceleration bias
    stdab = 0.5*[1, 1, 1]*1e-2;
    % Std turn-rate bias
    stdwb = 2*1.5* [1, 1, 1]*1e-3;

    stdab_mat = diag(stdab)^2;
    stdwb_mat = diag(stdwb)^2;

    if strcmp(str, 'lin')
        cov = dt*stdab_mat;
    elseif strcmp(str, 'ang')
        cov = dt*stdwb_mat;
    end
    bias_n = bias + [mvnrnd(0,cov(1,1)) mvnrnd(0,cov(2,2)) mvnrnd(0,cov(3,3))]';
end