 
%% Run the local filters individually

 for i = 1 : length(veh_ids)    

    veh_id = veh_ids(i);
    agent = agents{veh_id};

    % Assemble external and internal obs.
    H = [agent.H_int; agent.H_ext]; 
    z = [agent.z_int; agent.z_ext]; 
    R_add = diag( ones(length(agent.z_ext), 1)*system_setting.sigma_UWB^2 );
    R_int = agent.R_int;
    R = [R_int, zeros(size(R_int, 1), size(R_add, 2));  zeros(size(R_add, 1), size(R_int, 1)), R_add];
    h = afun_concat(agent.h_int, agent.h_ext);

    % Extended Kalman filter
    P = agent.P;
    %P = eye(size(P, 1));
    [x, P] = ekf_predict(agent.f, agent.F, agent.x, P, agent.Q);
    [x, P] = ekf_update(x, z, h, H, P, R);

    % Unscented Kalman filter
%     [x, P, X1, X2] = ukf_predict(f, x, P, Q);
%     [x, P]  = ukf_update(x, z, h, P, R, X1, X2);

    agents{veh_id}.apply_update(x, P);
end