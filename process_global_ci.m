%% Related materials:

% [1] Overview of federated filtering: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"
% [2] Original F-EKF paper (1988): https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=195473
% [3] Investigation on information sharing factor: https://pubs.acs.org/doi/pdf/10.1021/ie0511175
% [4] Proposed information sharing factor based on local filters covariances: https://pdfs.semanticscholar.org/5eea/21a9db2448c07234b99c94ed6015c27d643c.pdf
% [5] Information-sharing factor based on median and predicated states: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"

%% Run the local filters individually
for i = 1 : length(veh_ids)    

    veh_id = veh_ids(i);
    agent = agents{veh_id};
    if isempty(agent), continue; end

    % Assemble external and internal obs.
    H = [agent.H_int; agent.H_ext]; 
    z = [agent.z_int; agent.z_ext]; 
    R_add = diag( ones(length(agent.z_ext), 1)*system_setting.sigma_UWB^2 );
    R_int = agent.R_int;
    R = [R_int, zeros(size(R_int, 1), size(R_add, 2));  zeros(size(R_add, 1), size(R_int, 1)), R_add];
    h = afun_concat(agent.h_int, agent.h_ext);

    % Extended Kalman filter
     [x, P] = ekf_predict(agent.f, agent.F, agent.x, agent.P, agent.Q);
     [x, P] = ekf_update(x, z, h, H, P, R);

    % Unscented Kalman filter
%        [x, P, X1, X2] = ukf_predict(agent.f, agent.x, agent.P, agent.Q);
%        [x, P]  = ukf_update(x, z, h, P, R, X1, X2);

    agent.x = x;
    agent.P = P;
    %agent.apply_update(x, P);

end    

%% Covariance intersection

Ps = {};
Is = {};
for i = 1 : length(veh_ids)        
    Ps{i} = agents{veh_ids(i)}.P;
    Is{i} = inv(agents{veh_ids(i)}.P);
end

w0      = ones(length(veh_ids), 1) / length(veh_ids);
Aeq     = ones(1, length(veh_ids));
beq     = 1;
lb      = zeros(length(veh_ids), 1);
ub      = ones(length(veh_ids), 1);    
opts    = optimoptions('fmincon','Display','off');
w       = fmincon(@(x) trace(inv(ci_obj_func(x, Is))), w0, [], [], Aeq, beq, lb, ub, [], opts);   

%% Collect x and P, and fuse them
n   = size(Ps{1}, 1);
I_f = zeros(n, n);
x_f = zeros(n, 1);
I_avg = zeros(n, n);
for i = 1 : length(w)        
    I_i = Is{i};        
    I_f = I_f + w(i) * I_i;
    x_f = x_f + w(i) * I_i * agents{veh_ids(i)}.x;
    I_avg = I_avg + I_i;        

end
P_f = inv(I_f);
x_f = P_f * x_f; 

I_avg = I_avg / length(w);
P_avg = inv(I_avg);    

% Covariance visualization
if ~isempty(agents{sel_veh_id})
    figure(3); clf; hold on; 
    midxs = agents{sel_veh_id}.idxs;
    midxs = midxs(1:2);  
    links = agents{sel_veh_id}.links;
    for j = 1 : size(veh_ids, 1)               
        if min(eig(Ps{j})) < 0
            error(sprintf('P is not positive definite! Link: %i -> %i', links(j, 1), links(j, 2)));
        else
            error_ellipse(Ps{j}(midxs,midxs), [0 0], 'style', 'b');                        
        end      
    end
    error_ellipse(P_f(midxs,midxs), [0 0], 'style', 'r--');
    error_ellipse(P_avg(midxs,midxs), [0 0], 'style', 'g');
    axis equal;
end

% Filter reset
for i = 1 : length(veh_ids) 
    %agents{veh_ids(i)}.x = x_f;
    %agents{veh_ids(i)}.P = P_f;
    agents{veh_ids(i)}.apply_update(x_f, P_f);
end
