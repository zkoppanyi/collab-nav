clear variables;

% Clear persistent variables
clear CoopAgent

load('problem');

%% Related materials:

% [1] Overview of federated filtering: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"
% [2] Original F-EKF paper (1988): https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=195473
% [3] Investigation on information sharing factor: https://pubs.acs.org/doi/pdf/10.1021/ie0511175
% [4] Proposed information sharing factor based on local filters covariances: https://pdfs.semanticscholar.org/5eea/21a9db2448c07234b99c94ed6015c27d643c.pdf
% [5] Information-sharing factor based on median and predicated states: "Federated Filtering Revisited: New Directions to Distributed Systems Estimation and Filtering – a Case Study"
% [6] Original covariance intersection paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=609105



%% Settings 

system_setting.AOI = [750 1100; 1250 1600];                                 % Area of interest
sel_veh_id  = 20;                                                           % Selected vehicle ID
viz_frame_rate = 0.1;                                                       % Frame rate for visualization
dt = 0.05;

% Error characterization
simulation_settings;

%% Initialization
roads = problem.roads;
% plot(problem.vehicles(:, 2), problem.vehicles(:, 3), 'k.')
% axis equal;

% Get timestamps
ts = problem.vehicles(:,1);
ts = unique(ts);
form_t = 100; to_t = 160;
dxy_all = [];

% Discover how many vehicles we have in the give time interval
veh_ids_total = [];
for t = form_t : dt*5 : to_t
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);    
    veh_ids = unique(epoch(:,2));
    veh_ids_total = [veh_ids_total; veh_ids]; 
end
veh_ids_total = unique(veh_ids_total);
agents{length(veh_ids_total)} = {};

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Algorithm

% Init filter
n_veh = length(veh_ids_total);
n_states = 4;
x_m = zeros(n_veh*n_states, 1);
P_m = eye(n_veh*n_states, n_veh*n_states);

master_agent = CoopAgent(-1, [0 0 0 0], eye(4), 0, n_veh, system_setting);

init_rewind = form_t;
init_i = 0;
init_time = 0.1/dt;
is_init = 0;
weights = [];

for t = form_t : dt : to_t

    if init_i < init_time
        t = init_rewind;
        init_i = init_i + 1;
        is_init = 1;
    else
        is_init = 0;
    end
    
    fprintf("t= %.1f\n", t);     
    
    idx = find( abs(problem.vehicles(:, 1) - t) < dt/2 );
    epoch = problem.vehicles(idx, :);
    
    if isempty(epoch)
        continue;
    end
    
    % Find trajectory within the AOI
    epoch = epoch(and(and(and(epoch(:,3) > system_setting.AOI (1,1), epoch(:,3) < system_setting.AOI (1,2)), epoch(:,4) > system_setting.AOI (2,1)), epoch(:,4) < system_setting.AOI (2,2)), :);   
           
    % Vehicles at the current epoch
    veh_ids = unique(epoch(:,2));       
    if isempty(find(veh_ids == sel_veh_id, 1))
        continue;
    end
    
    % New cycle; reset internal variables
     for i = 1 : length(agents)
        if isempty(agents{i}), continue; end            
        agents{i}.reset();
    end
        
    %% Internal observations    
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);      
                
        % Get internal observations
        inter_obs = epoch(epoch(:,2) == veh_id, :);
                 
        % Get internal observation
        xy = inter_obs(3:4)';
        b  = inter_obs(5)/180*pi;
        v  = inter_obs(6);
        upt = [xy', v, b];  
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = CoopAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
           
            % Initialized master filter x and P
            %idxs = agents{veh_id}.idxs;
            %master_agent.x(idxs) = agents{veh_id}.x(idxs);
            %master_agent.P(idxs, idxs) = agents{veh_id}.P(idxs, idxs);   
            master_agent.comm(agents{veh_id} , 'share-states');  
            init_i = 0; % apply intialzation
            init_rewind = t;
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();         
       
        % Measurement\observation        
        if is_init == 1
%            agents{veh_id}.build_int_update('GPS', upt); % init mode
            agents{veh_id}.build_int_update('GPS+IMU', upt);
        else       
            is_gps = ~isempty(find(gps_agent == veh_id));
            if and(mod(t, 0.2) == 0, is_gps)
                agents{veh_id}.build_int_update('GPS+IMU', upt);
            else
                agents{veh_id}.build_int_update('IMU', upt);            
            end
        end
                
    end    
    
    %% External observations: relative ranging       
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);        
        xy_gt = agents{veh_id}.gt(end, :);
        
        fprintf('%i -> ', veh_id);
        
        % Add relative ranging
        nidx = find( and( sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) < system_setting.com_radius, epoch(:, 2) ~= veh_id )  ); 
        links_loc = epoch(nidx, :);                                         % neighbors
        
        for j = 1 : size(links_loc, 1)
            
            % Get indeces of neighbor
            veh_id_j = links_loc(j, 2);
            
            % This means that the vehicle is out of AOI
            if isempty(agents{veh_id_j})
                break;
            end                      
            
            agents{veh_id}.comm(agents{veh_id_j}, 'share-states');
            agents{veh_id}.comm(agents{veh_id_j}, 'UWB');            
            fprintf('%i ', veh_id_j);                        
        end        
        fprintf('\n');
        
        % Add infrastructure ranging
        if ~isempty(infra_nodes)
            nidx = find( sqrt( (infra_nodes(:, 2) - xy_gt(1)).^2 + (infra_nodes(:, 3) - xy_gt(2)).^2 ) < system_setting.com_radius ); 
            for j = 1 : length(nidx)
                agents{veh_id}.comm(infra_nodes(nidx(j), :), 'v2i');
            end
        end
        
    end
    fprintf('\n');
 
    %% Federated fitlering: calculate beta
    
    % Calculate DOP
%      A = [];
%     for i = 1 : length(agents)        
%         if isempty(agents{i}), continue; end
%          
%         agent = agents{i};
%         links = agent.links;
%         for j = 1 : size(links, 1)
%             if (links(j, 2) < 1000)
%             
%                 agent2 = agents{links(j, 2)};           
%                 sidx_i = agent.sidx;
%                 sidx_j = agent2.sidx;
% 
%                 A_line = zeros(1, agent.max_neighbors * agent.n_states); 
%                 x1 = agent.x(sidx_i);       y1 = agent.x(sidx_i+1); 
%                 x2 = agent2.x(sidx_j);      y2 = agent2.x(sidx_j+1); 
%                 r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
% 
%                 A_line(sidx_i)        = (x1-x2)/r;
%                 A_line(sidx_i+1)      = (y1-y2)/r;
%                 A_line(sidx_j)        = (x2-x1)/r;
%                 A_line(sidx_j+1)      = (y2-y1)/r;
%             
%             else
%                 
%                 sidx_i = agent.sidx;
%                 iidx = infra_nodes(:,1) == links(j, 2);
%                 infra_xy = infra_nodes(iidx, 2:3);
%     
%                 A_line = zeros(1, agent.max_neighbors * agent.n_states); 
%                 x1 = agent.x(sidx_i);  y1 = agent.x(sidx_i+1); 
%                 x2 = infra_xy(1);      y2 = infra_xy(2); 
%                 r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
% 
%                 A_line(sidx_i)        = (x1-x2)/r;
%                 A_line(sidx_i+1)      = (y1-y2)/r;
%                 
%             end
%             A = [A; A_line];    
%         end           
%         
%         
%     end      
%     DOP = pinv(A'*A);
    
    % Calcualte weights
%     weights = []; 
%     for i = 1 : length(agents)        
%         if isempty(agents{i}), continue; end
%         idxs = agents{i}.idxs;
%         DOP_l = DOP(idxs, idxs);
%         w = sqrt(DOP_l(1,1)^2 + DOP_l(2,2)^2);
%         nn = size(agents{i}.links, 1);
%         weights = [weights; i, w];
%     end  
%     
%     
%     %betas3 = ones(25, 1);
%     if ~isempty(weights)
%         %betas = weights(:, 2);
%         %betas = weights(:, 2) ./ (sum(weights(:, 2))) * 1e-3;
%         betas = ones(size(weights, 1), 1)*1e-6;
%         %betas = weights(:, 2) ./ (sum(weights(:, 2)) + length(weights));
%         %betas = ones(size(traces,1), 1) / size(traces, 1);
%         %betas = ones(size(traces,1), 1);
%         for i = 1 : size(weights, 1) 
% 
%             x = agents{weights(i, 1)}.x;
%             P = agents{weights(i, 1)}.P;
%             Q = agents{weights(i, 1)}.Q;
% 
%             b2 = betas(i);
%             if betas(i) ~= 0
%                 %P = (1/betas(i))*P * ((length(betas) + sum(weights(:, 2))) * P);
%                 %Q = (1/betas(i))*Q * ((length(betas) + sum(weights(:, 2))) * Q);
%                 
%                 i2 = i + size(weights, 1) / 2;
%                 %P = (1/betas(i2))*P*((1/betas(i))*P);
%                 P = (1/betas(i))*P;
%                 %P = 1/betas(i+1)*P;
%                 %betas3(weights(i, 1)) = betas3(weights(i, 1)) * (1/betas(i));
%                 %Q = (1/betas(i))*Q;
%                 
%             end
% 
%             % Update agents' states       
%             agents{weights(i, 1)}.P = P;       
%             %agents{weights(i, 1)}.Q = Q;
% 
%         end 
%     end
    
    % Covariance intersection
    Ps = {};
    for i = 1 : length(veh_ids)        
        agent = agents{veh_ids(i)};
        Ps{i} = agent.P;
    end   
    
    %ci_obj_func(ones(length(veh_ids), 1), Ps)
    x0 = ones(length(veh_ids), 1) / length(veh_ids) + 0.2;
    %x0 =  [0.1109    0.0000    0.0000    0.0000    0.0001    0.0000    0.0000    0.0001    0.0000];
    %x0 = betas;
    Aeq = ones(length(veh_ids));
    beq = ones(length(veh_ids), 1);
    lb = zeros(length(veh_ids), 1);
    ub = ones(length(veh_ids), 1);    
    betas2 = fmincon(@(x) trace(ci_obj_func(x, Ps)), x0, [], [], Aeq, beq, lb, ub);
    if is_init 
        betas2 = betas2 * 1e-3
    end
    
    for i = 1 : length(veh_ids)        
        agents{veh_ids(i)}.P = (1/betas2(i))*agent.P;
        %agents{veh_ids(i)}.Q = (1/betas2(i))*agent.Q;
    end
    
    %return
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
    
    master_agent.build_predict();

    [x_m, P_m] = ekf_predict(master_agent.f, master_agent.F, master_agent.x, master_agent.P, master_agent.Q);    
    master_agent.x = x_m;
    master_agent.P = P_m;
 
     
    %% Collect x and P, and fuse them
    P_f = inv(master_agent.P);
    x_f = inv(master_agent.P) * master_agent.x;
    for i = 1 : length(agents)
        if isempty(agents{i}), continue; end
        x = agents{i}.x;
        P = agents{i}.P;
        
        P_inv = inv(P);
        P_f = P_f + P_inv;
        x_f = x_f + P_inv * x;
    end
    P_f = inv(P_f);
    x_f = P_f * x_f;
    master_agent.apply_update(x_f, P_f);
    
    % Filter reset
    for i = 1 : length(agents)
        if isempty(agents{i}), continue; end
        agents{i}.apply_update(x_f, P_f);
    end
    
    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
    if is_init
        title('Initialization');
    end
    plot(roads(:, 1), roads(:, 2), 'k.');  
    if ~isempty(infra_nodes)
        plot(infra_nodes(:, 2), infra_nodes(:, 3), 'r.', 'MarkerSize', 15);    
    end
    for i = 1 : length(agents)        
        if isempty(agents{i}), continue; end
        
        x_hist    = agents{i}.my_x_hist();
        x_loc     = x_hist(end, :);
        xy_loc_gt = agents{i}.gt(end, :);
        z_loc     = agents{i}.z_int_hist(end, :);
        z_loc     = z_loc(agents{i}.sidx:(agents{i}.sidx+1));

        plot(x_loc(1), x_loc(2), 'b.', 'MarkerSize', 15);
        plot(xy_loc_gt(1), xy_loc_gt(2), 'g.', 'MarkerSize', 15);  
        plot(z_loc(1), z_loc(2), 'r.', 'MarkerSize', 15); 

        if agents{i}.id == sel_veh_id
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w', 'Color', 'r');
        else
            text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w');
        end
        
        if ~isempty(agents{i}.links)
            links_agents = agents{i}.links(agents{i}.links(:,2) < 1000, :);
            for j = 1 : size(links_agents, 1)
                nidx = links_agents(j, 2);
                xy_loc_gt_j = agents{nidx}.gt(end, :);
                plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'g-');
            end

            links_infra = agents{i}.links(agents{i}.links(:,2) >= 1000, :);
            for j = 1 : size(links_infra, 1)
                nidx = find( links_infra(j, 2) == infra_nodes(:,1) );
                xy_loc_gt_j = infra_nodes(nidx, 2:3);
                plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'r.-');
            end
       end
    end

    grid on; axis equal;    
    xlim(system_setting.AOI(1,1:2)); ylim(system_setting.AOI(2,1:2));       
    set(gca, 'FontSize', 12);    
    xlabel('[m]'); ylabel('[m]');
    
    % Selected vehicle's internal states    
    agent = agents{sel_veh_id};
    if ~isempty(agent)
        
        x_hist = agent.my_x_hist();
        
        idxs = agent.idxs;
        x_hist_master = master_agent.x_hist(:, idxs);
        
        n = size(x_hist, 1);
        dxy = (x_hist_master(end,:)-agent.gt(end, :))';

        figure(2); clf; hold on;                

        subplot(4, 1, 1); hold on; 
        plot(1:n, x_hist(:,1), 'b.-' );
        plot(1:n, agent.gt(:, 1), 'g.-' );
        plot(1:n,  x_hist_master(:,1), 'r.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 2); hold on; 
        plot(1:n,x_hist(:,2), 'b.-' );
        plot(1:n, agent.gt(:, 2), 'g.-' );
        plot(1:n,  x_hist_master(:,2), 'r.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 3); hold on; 
        plot(1:n, x_hist(:,3), 'b.-' );
        plot(1:n, agent.gt(:, 3), 'g.-' );
        plot(1:n,  x_hist_master(:,3), 'r.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 4); hold on; 
        plot(1:n, x_hist(:,4)/pi*180, 'b.-' );
        plot(1:n, agent.gt(:, 4)/pi*180, 'g.-' );
        plot(1:n,  x_hist_master(:,4)/pi*180, 'r.-' );
        title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
        set(gca, 'FontSize', 10);       
    end
        
    
    % Show sparisty of P
%     figure(3);
%     show_sparisty_pattern(P, 3);
    
    % Compate trajectries to ground truth   
    fprintf('Time: %.1f\n', t);    
    dxy_all = [];
    for i = 1 : length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            x_hist =  agents{veh_ids(i)}.my_x_hist();
            %dxy = agents{veh_ids(i)}.gt(:, 1:2) - x_hist(:, 1:2);
            dxy = agents{veh_ids(i)}.gt(end, 1:2) - x_hist(end, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            %dxy_all  = [dxy_all; dxy(end)];
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL MEAN: [%.3f] MEDIAN: [%.3f]\n', mean(dxy_all), median(dxy_all));
    
    pause(viz_frame_rate);
end



