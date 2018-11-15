if ~exist('is_batch_process', 'var')
    clear variables;   
    method = 4;
end

% Clear persistent variables
clear CoopAgent
clear CoopAgent2
clear agents;

% Load problem
load('problem');

simulation_settings;
    
%% Methods' settings
% 1 - Tracking
% 2 - Naive state sharing
% 3 - Global CI
% 4 - Local CI
% 5 - Consensus

if method == 1
    
    m_name = 'tracking';
    is_share_states = 0;
    init_time = 0.1/dt;
    
elseif method == 2
    
    m_name = 'state_sharing';
    is_share_states = 1;   
    init_time = 0.1/dt;
    
elseif method == 3
    
    m_name = 'global_ci';
    is_share_states = 0;
    init_time = 0.2/dt;
    
elseif method == 4
    
    m_name = 'local_ci';
    is_share_states = 0;
    init_time = 0.5/dt;
    
elseif method == 5
    
    m_name = 'consensus';
    is_share_states = 0;
    init_time = 0.5/dt;
    
end

    
%% Initialization

roads = problem.roads;

% Get timestamps
ts = problem.vehicles(:,1);
ts = unique(ts);
form_t = 100; to_t = 110;
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
veh_ids_total = [10 16 8 20 13 6 17];  % fixed ids
%agents{length(veh_ids_total)} = {};
agents{250} = {};

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Simulation Algorithm

% Init filter
n_veh = length(veh_ids_total);
n_states = 4;
x_m = zeros(n_veh*n_states, 1);
P_m = eye(n_veh*n_states, n_veh*n_states);

% For initialization periods
init_rewind = form_t;
init_i = 0;
is_init = 0;
t = form_t;
iter_i = 0;

prev_n_veh_ids = 0;
prev_n_groups = 0;
while t <= to_t

    if init_i < init_time
        t = init_rewind;
        init_i = init_i + 1;
        is_init = 1;
    else
        is_init = 0;
        t = t + dt;
        iter_i = iter_i + 1;
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
    if prev_n_veh_ids ~= length(veh_ids)
        % A vehicle left or joined: apply intialzation
        %is_init = 1;
        %init_i = 0; 
        %init_rewind = t;
    end
    prev_n_veh_ids = length(veh_ids);
    
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end
    
    % Fixed number of vehicles
    rm_idx = [];
    rm_epoch = [];
    for i = 1 : length(veh_ids)
        if isempty(find(veh_ids(i) == veh_ids_total))
            rm_idx = [rm_idx; i]; 
            idx = find(veh_ids(i) == epoch(:, 2));
            rm_epoch = [rm_epoch; idx];
        end
    end
    veh_ids(rm_idx) = [];
    epoch(rm_epoch, :) = [];
    
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
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = CoopAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
           %agents{veh_id} = CoopAgent2(veh_id, x_init, P_init, t, n_veh, system_setting);
           
           % A new vehicle joins: apply intialzation
           is_init = 1;
           init_i = 0; 
           init_rewind = t;
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();         
       
        % Inter-node measurement\observation update
        upt = [xy', v, b, t, is_init];  

        if is_init == 1
%            agents{veh_id}.build_int_update('GPS', upt); % init mode
            agents{veh_id}.build_int_update('ZUPT', upt);
        else       
            is_gps = ~isempty(find(gps_agent == veh_id));
            if and(mod(iter_i, 4) == 0, is_gps)
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
        if simulation_scenario == 4
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
                
            % State sharing if applicable: only local internal states shared
            if is_share_states == 1
                agents{veh_id}.comm(agents{veh_id_j}, 'share-states');
                agents{veh_id_j}.comm(agents{veh_id}, 'share-states');
            end

            % Apply UWB update only, if this vehilce didn't do earlier
            if ~isempty(agents{veh_id}.links)
                if isempty(find( agents{veh_id}.links(:, 2) == veh_id_j ))
                    agents{veh_id}.comm(agents{veh_id_j}, 'UWB'); 
                end
            else
                agents{veh_id}.comm(agents{veh_id_j}, 'UWB'); 
            end
            
            % Apply UWB update only, if other vehilce didn't do so
            if ~isempty(agents{veh_id_j}.links)
                if isempty(find( agents{veh_id_j}.links(:, 2) == veh_id ))
                    agents{veh_id_j}.comm(agents{veh_id}, 'UWB');            
                end
            else
                agents{veh_id_j}.comm(agents{veh_id}, 'UWB');
            end
            
            fprintf('%i ', veh_id_j);  
            
        end        
        fprintf('\n');
        
        % Add infrastructure ranging
        if simulation_scenario ~= 4
            if ~isempty(infra_nodes)
                nidx = find( sqrt( (infra_nodes(:, 2) - xy_gt(1)).^2 + (infra_nodes(:, 3) - xy_gt(2)).^2 ) < system_setting.com_radius ); 
                for j = 1 : length(nidx)
                    agents{veh_id}.comm(infra_nodes(nidx(j), :), 'v2i');
                end
            end
        end
        
    end       
    fprintf('\n');
  
     if simulation_scenario == 4
         for i = 1 : size(infra_nodes, 1)
             dist = sqrt( (infra_nodes(i, 2) - epoch(:, 3)).^2 + (infra_nodes(i, 3) - epoch(:, 4)).^2 );
             [~, nidx] = min(dist);
                agents{epoch(nidx, 2)}.comm(infra_nodes(i, :), 'v2i');
         end
     end
     
     % Make the communication graph fully connected
%      [~, G, ~] = get_local_conn(agents, veh_ids, veh_ids(1));
%      gr = conncomp(G);
%      ugr = unique(gr);
%      if length(ugr) > 1
%          ref_gr = mode(gr);         
%          idx = find(ref_gr == gr);
%          ref_node = agents{veh_ids(idx(1))};
%          for k = 1 : length(ugr)
%              if ugr(k) ~= ref_gr
%                 idx = (ugr(k) == gr);                
%              end
%          end
%      end
     
%       [~, G, ~] = get_local_conn(agents, veh_ids, veh_ids(1));
%       gr = conncomp(G);
%       ugr = unique(gr);
%       if prev_n_groups > length(ugr)
%           
%            % Group joins: apply intialization
%            is_init = 1;
%            init_i = 0; 
%            init_rewind = t;           
%       end
%       prev_n_groups = length(ugr);
    
%      % Initialize all vehciles: faster intialization
%      % this step can be implemented as gossip algorithm
%      if is_init
%         for i = 1 : length(veh_ids) 
%             for j = 1 : length(veh_ids) 
%                 agents{veh_ids(i)}.comm(agents{veh_ids(j)}, 'share-states');
%             end
%         end
%      end
  
   %% Methods
   
   if or(method == 1, method == 2)
        process_tracking;
   elseif method == 3
        process_global_ci;
    elseif method == 4
        process_local_ci;        
   elseif method == 5
        process_consensus;       
   end

    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
     if is_init
        title(sprintf("init t= %.3f\n", t));
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
    if is_init == 0
    
        agent = agents{sel_veh_id};
    
        if ~isempty(agent)

            x_hist = agent.my_x_hist();
            gt = agent.gt;

            % remove init periods
            idx = find(gt(:, end) == 0);
            x_hist = x_hist(idx, :);
            gt = gt(idx, :);

            n = size(x_hist, 1);
            dxy = (x_hist(end,1:agent.n_states)-agent.gt(end, 1:agent.n_states))';

            figure(2); clf; hold on;                

            subplot(4, 1, 1); hold on; 
            plot(1:n, x_hist(:,1), 'b.-' );
            plot(1:n, gt(:, 1), 'g.-' );
            title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
            set(gca, 'FontSize', 10);

            subplot(4, 1, 2); hold on; 
            plot(1:n,x_hist(:,2), 'b.-' );
            plot(1:n, gt(:, 2), 'g.-' );
            title(sprintf('Y= %.3f', dxy(2)));
            set(gca, 'FontSize', 10);

            if agent.n_states >= 4
                subplot(4, 1, 3); hold on; 
                plot(1:n, x_hist(:,3), 'b.-' );
                plot(1:n, gt(:, 3), 'g.-' );
                title(sprintf('v= %.3f', dxy(3)));
                set(gca, 'FontSize', 10);

                subplot(4, 1, 4); hold on; 
                plot(1:n, x_hist(:,4)/pi*180, 'b.-' );
                plot(1:n, gt(:, 4)/pi*180, 'g.-' );
                title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
                set(gca, 'FontSize', 10);       
            end
        end
    end
        
    
    % Show sparisty of P
%     figure(3);
%     show_sparisty_pattern(P, 3);
    
    % Compate trajectries to ground truth   
    fprintf('Time: %.1f Method: %s\n', t, m_name);    
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
    
    figure(1); 
     if ~is_init
        title(sprintf("t= %.3f err=%.3f\n", t, mean(dxy_all)));
     end
    
    pause(viz_frame_rate);
end


sol.agents = agents;
sol.system_setting = system_setting;
sol.problem = problem;
save(['results\' sprintf('solution_%s_%i', m_name, simulation_scenario)], 'sol')

