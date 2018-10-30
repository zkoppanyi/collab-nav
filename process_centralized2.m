clear variables;

%% Settings 
simulation_settings;
 
%gps_error_tests = [0.1 0.5 1.0 5.0 10.0];
gps_error_tests = 10.0;

for test_i = 1 : length(gps_error_tests)
    system_setting.sigma_GPS = gps_error_tests(test_i);
    
    load('problem');
    clear agents

    % Clear persistent variables
    clear SelfishAgent
    clear SelfishAgent2
    
%% Initialization
roads = problem.roads;
% plot(problem.vehicles(:, 2), problem.vehicles(:, 3), 'k.')
% axis equal;

% Get timestamps
ts = problem.vehicles(:,1);
ts = unique(ts);

form_t = 100; to_t = 130;

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
%agents{length(veh_ids_total)} = {};
agents{250} = {};

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Algorithm

% Init filter
n_veh = length(veh_ids_total);
n_states = 2;
x = zeros(n_veh*n_states, 1);
P = eye(n_veh*n_states, n_veh*n_states);

for t = form_t : dt : to_t

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
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end


     if isempty(veh_ids)
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
           %x_init = [xy; v; b];
           %P_init = eye(4,4);  
           x_init = xy;
           P_init = eye(2);  
           agents{veh_id} = SelfishAgent2(veh_id, x_init, P_init, t, n_veh, system_setting);
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();     
     
        % Measurement\observation       
        is_gps = ~isempty(find(gps_agent == veh_id));
        if and(mod(t, 0.2) == 0, is_gps)
            agents{veh_id}.build_int_update('GPS+IMU', upt);
        else
            agents{veh_id}.build_int_update('IMU', upt);            
        end
                
    end    
    
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);        
        xy_gt = agents{veh_id}.gt(end, :);
        
        fprintf('%i -> ', veh_id);
        
        % Add relative ranging
        if simulation_scenario == 3
            dist = sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) ; 
            [vals, idx] = sort(dist);
            links_loc = epoch(idx(2:4), :);                                    % neighbors            
        else
            nidx = find( and( sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) < system_setting.com_radius, epoch(:, 2) ~= veh_id )  ); 
            links_loc = epoch(nidx, :);                                         % neighbors        
        end

        for j = 1 : size(links_loc, 1)
            
            % Get indeces of neighbor
            veh_id_j = links_loc(j, 2);
            
            % This means that the vehicle is out of AOI
            if isempty(agents{veh_id_j})
                continue;
            end                      
            
            %agents{veh_id}.comm(agents{veh_id_j}, 'share-states');
            agents{veh_id}.comm(agents{veh_id_j}, 'UWB');            
            fprintf('%i ', veh_id_j);                        
        end        
        fprintf('\n');
        
        % Add infrastructure ranging
        if simulation_scenario ~= 3
            if ~isempty(infra_nodes)
                nidx = find( sqrt( (infra_nodes(:, 2) - xy_gt(1)).^2 + (infra_nodes(:, 3) - xy_gt(2)).^2 ) < system_setting.com_radius ); 
                for j = 1 : length(nidx)
                    agents{veh_id}.comm(infra_nodes(nidx(j), :), 'v2i');
                end
            end
        end
        
    end       
    fprintf('\n');
  
     if simulation_scenario == 3
         for i = 1 : size(infra_nodes, 1)
             dist = sqrt( (infra_nodes(i, 2) - epoch(:, 3)).^2 + (infra_nodes(i, 3) - epoch(:, 4)).^2 );
             [~, nidx] = min(dist);
             agents{epoch(nidx, 2)}.comm(infra_nodes(i, :), 'v2i');
         end
     end
    %% Collect data (agents transmitting states to master)
    
    % Init filter matrices
    Q = eye(n_veh*n_states, n_veh*n_states);
    R = eye(n_veh*n_states, n_veh*n_states);
    F = eye(n_veh*n_states, n_veh*n_states);
    f = afun_create(n_veh*n_states, @(x) 0);
    H_int = zeros(n_veh*n_states, n_veh*n_states);
    h_int = afun_create(n_veh*n_states, @(x) 0);
    z_int = zeros(n_veh*n_states, 1);
    
    % Internal observations        
    for i = 1 : length(agents)

        agent = agents{i};
        if isempty(agent), continue; end

        % Indeces
        sidx = (agent.sid-1)*n_states+1;
        idxs = sidx:(sidx+n_states-1);

        % Collect data and form matrices
        
        x(idxs)         = agent.x;
        P(idxs, idxs)   = agent.P; 
        F(idxs, idxs)   = agent.F;
        f{idxs(1)}      = @(x) agent.f{1}(x(idxs));
        f{idxs(2)}      = @(x) agent.f{2}(x(idxs));
        %f{idxs(3)}      = @(x) agent.f{3}(x(idxs));
        %f{idxs(4)}      = @(x) agent.f{4}(x(idxs));
        Q(idxs, idxs)   = agent.Q;
        R(idxs, idxs)   = agent.R_int;                
     
        % Measurement\observation matrix            
        H_int(idxs, idxs) = agent.H_int;
        z_int(idxs)       = agent.z_int;
        h_int{sidx}       = @(x) agent.h_int{1}(x(idxs));
        h_int{sidx+1}     = @(x) agent.h_int{2}(x(idxs));
        %h_int{sidx+2}     = @(x) agent.h_int{3}(x(idxs));
        %h_int{sidx+3}     = @(x) agent.h_int{4}(x(idxs));   
        
    end
    
    % External observations    
    H_ext = [];
    h_ext = {};
    z_ext = [];    
    for i = 1 : length(agents)
        
        agent = agents{i};
        if isempty(agent), continue; end

        H_line = agent.H_ext;
        H_ext = [H_ext; H_line];           %#ok
        z_ext = [z_ext; agent.z_ext];      %#ok            
        h_ext = afun_concat( h_ext, agent.h_ext);                
    end
    
    % Let's assemble the problem
    H = [H_int; H_ext]; 
    z = [z_int; z_ext]; 
    R_add = diag( ones(length(z_ext), 1)*system_setting.sigma_UWB^2 );
    R_int = R;
    R = [R, zeros(size(R, 1), size(R_add, 2));  zeros(size(R_add, 1), size(R, 1)), R_add];  %#ok
    h = afun_concat(h_int, h_ext);
    
    %% Centralized filtering
    
    % Extended Kalman filter        
    [x, P] = ekf_predict(f, F, x, P, Q); 
    [x, P] = ekf_update(x, z, h, H, P, R);
    
    % Unscented Kalman filter
%     [x, P, X1, X2] = ukf_predict(f, x, P, Q);
%     [x, P]  = ukf_update(x, z, h, P, R, X1, X2);


    %% Transmit back the estimates    
   for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);           
        agent = agents{veh_id};

        % Indeces
        sidx = (agent.sid-1)*n_states+1;
        idxs = sidx:(sidx+n_states-1);

        % Update agents' states        
        agent.apply_update(x(idxs), P(idxs, idxs));
        
   end     
    
   
    % Covariance visualization
    figure(3); clf; hold on; 
    Ps = agents{sel_veh_id}.P;
    for j = 1 : size(veh_ids, 1)               
        if min(eig(P)) < 0
            error(sprintf('P is not positive definite! Link: %i -> %i', links(j, 1), links(j, 2)));
        else
            error_ellipse(Ps(1:2,1:2), [0 0], 'style', 'b');                        
        end      
    end
    axis equal;
    
    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
    %title(sprintf('t= %.1f Acc: %.2f', t, system_setting.sigma_GPS));
    plot(roads(:, 1), roads(:, 2), 'k.');   
    if ~isempty(infra_nodes)
        plot(infra_nodes(:, 2), infra_nodes(:, 3), 'r.', 'MarkerSize', 25);    
    end
    for i = 1 : length(agents)        
        if isempty(agents{i}), continue; end
        
        x_loc     = agents{i}.x;
        xy_loc_gt = agents{i}.gt(end, :);
        z_loc     = agents{i}.z_int_hist(end, :);

       

%         if agents{i}.id == sel_veh_id
%             text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w', 'Color', 'r');
%         else
%             text(xy_loc_gt(1), xy_loc_gt(2)+20, sprintf('%i',  agents{i}.id), 'BackgroundColor', 'w');
%         end
        
        if ~isempty(agents{i}.links)
            links_agents = agents{i}.links(agents{i}.links(:,2) < 1000, :);
%             for j = 1 : size(links_agents, 1)
%                 nidx = links_agents(j, 2);
%                 xy_loc_gt_j = agents{nidx}.gt(end, :);
%                 plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'g-', 'LineWidth', 1.5);
%             end

            links_infra = agents{i}.links(agents{i}.links(:,2) >= 1000, :);
            for j = 1 : size(links_infra, 1)
                nidx = find( links_infra(j, 2) == infra_nodes(:,1) );
                xy_loc_gt_j = infra_nodes(nidx, 2:3);
                plot([xy_loc_gt(1), xy_loc_gt_j(1)], [xy_loc_gt(2), xy_loc_gt_j(2)], 'r.-', 'LineWidth', 1.5);
            end
        end
        plot(x_loc(1), x_loc(2), 'b.', 'MarkerSize', 25);
        plot(xy_loc_gt(1), xy_loc_gt(2), 'g.', 'MarkerSize', 25);  
        %plot(z_loc(1), z_loc(2), 'r.', 'MarkerSize', 25); 
    end

    grid on; axis equal;    
    xlim(system_setting.AOI(1,1:2)); ylim(system_setting.AOI(2,1:2));       
    set(gca, 'FontSize', 14);    
    xlabel('[m]'); ylabel('[m]');
    
    
    % Selected vehicle's internal states    
    agent = agents{sel_veh_id};
    if ~isempty(agent)
        n = size(agent.x_hist, 1);
        dxy = (agent.x_hist(end, :)-agent.gt(end, :))';

        figure(2); clf; hold on;                

        subplot(4, 1, 1); hold on; 
        plot(1:n, agent.x_hist(:, 1), 'b.-' );
        plot(1:n, agent.gt(:, 1), 'g.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 2); hold on; 
        plot(1:n, agent.x_hist(:, 2), 'b.-' );
        plot(1:n, agent.gt(:, 2), 'g.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 3); hold on; 
        plot(1:n, agent.x_hist(:, 3), 'b.-' );
        plot(1:n, agent.gt(:, 3), 'g.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 4); hold on; 
        plot(1:n, agent.x_hist(:, 4)/pi*180, 'b.-' );
        plot(1:n, agent.gt(:, 4)/pi*180, 'g.-' );
        title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
        set(gca, 'FontSize', 10);       
    end
        
    
    % Show sparisty of P
%     figure(3);
%     show_sparisty_pattern(P, 3);
    
    % Compate trajectries to ground truth   
    fprintf('Time: %.1f\n', t);    
    for i = 1 : length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            dxy = agents{veh_ids(i)}.gt(:, 1:2) - agents{veh_ids(i)}.x_hist(:, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL MEAN: [%.3f] MEDIAN: [%.3f]\n', mean(dxy_all), median(dxy_all));
    
    pause(viz_frame_rate);
    %return
end

solution_centralized.agents = agents;
solution_centralized.system_setting = system_setting;
solution_centralized.problem = problem;
%save(sprintf('solution_centralized_%i_%i', simulation_scenario, gps_error_tests(test_i)*10), 'solution_centralized')
save(sprintf('solution_centralized_%i', simulation_scenario), 'solution_centralized')

end
