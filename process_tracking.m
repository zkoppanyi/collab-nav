clear variables;

% Clear persistent variables

%% Related materials:



%% Settings 

system_setting.AOI = [750 1100; 1250 1600];                                 % Area of interest
sel_veh_id  = 20;                                                           % Selected vehicle ID
viz_frame_rate = 0.1;                                                       % Frame rate for visualization
dt = 0.05;

% Error characterization
simulation_settings;

%gps_error_tests = [0.1 0.5 1.0 5.0 10.0];
gps_error_tests = [5.0];

for test_i = 1 : length(gps_error_tests)
    system_setting.sigma_GPS = gps_error_tests(test_i);
    
    load('problem');
    clear agents
    clear CoopAgent
    clear CoopAgent2
    
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
%veh_ids_total = [10 16 8 20 13 6 17];
%agents{length(veh_ids_total)} = {};
agents{250} = {};

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Algorithm

% Init filter
n_veh = length(veh_ids_total);
n_states = 4;
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
    
    % Fixed number of vehicles
    rm_idx = [];
    for i = 1 : length(veh_ids)
        if isempty(find(veh_ids(i) == veh_ids_total))
            rm_idx = [rm_idx; i];
        end
    end
    veh_ids(rm_idx) = [];
    
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end


    
    
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
        %upt = [xy'];  
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = CoopAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
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
    
    %% External observations: relative ranging
       
    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);        
        xy_gt = agents{veh_id}.gt(end, :);
        
        fprintf('%i -> ', veh_id);
        
        % Add relative ranging
        if simulation_scenario == 3
            dist = sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) ; 
            [vals, idx] = sort(dist);
            links_loc = epoch(idx(2:4), :);                                % neighbors            
        else
            nidx = find( and( sqrt( (epoch(:, 3) - xy_gt(1)).^2 + (epoch(:, 4) - xy_gt(2)).^2 ) < system_setting.com_radius, epoch(:, 2) ~= veh_id )  ); 
            links_loc = epoch(nidx, :);                                    % neighbors        
        end

        
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
    
    % Add infrastructure ranging
    if simulation_scenario == 3
         for i = 1 : size(infra_nodes, 1)
             dist = sqrt( (infra_nodes(i, 2) - epoch(:, 3)).^2 + (infra_nodes(i, 3) - epoch(:, 4)).^2 );
             [~, nidx] = min(dist);
             agents{epoch(nidx, 2)}.comm(infra_nodes(i, :), 'v2i');
         end
    end
     
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
        %P = agent.P;
        P = eye(size(agent.P, 1))*1;
        [x, P] = ekf_predict(agent.f, agent.F, agent.x, P, agent.Q);
        [x, P] = ekf_update(x, z, h, H, P, R);

        % Unscented Kalman filter
    %     [x, P, X1, X2] = ukf_predict(f, x, P, Q);
    %     [x, P]  = ukf_update(x, z, h, P, R, X1, X2);
    
        agents{veh_id}.apply_update(x, P);
    end

    
    %% Visualization and analytics
    
    % Vehicles positions and communication links
    figure(1); clf; hold on;
    title(sprintf('t= %.1f Acc: %.2f', t, system_setting.sigma_GPS));
    plot(roads(:, 1), roads(:, 2), 'k.');  
    if ~isempty(infra_nodes)
        plot(infra_nodes(:, 2), infra_nodes(:, 3), 'r.', 'MarkerSize', 15);
    end
    for i = 1 : length(agents)        
        if isempty(agents{i}), continue; end
        
        x_hist     = agents{i}.my_x_hist();
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
        n = size(x_hist, 1);
        dxy = (x_hist(end,:)-agent.gt(end, :))';

        figure(2); clf; hold on;                

        subplot(4, 1, 1); hold on; 
        plot(1:n, x_hist(:,1), 'b.-' );
        plot(1:n, agent.gt(:, 1), 'g.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 2); hold on; 
        plot(1:n,x_hist(:,2), 'b.-' );
        plot(1:n, agent.gt(:, 2), 'g.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 3); hold on; 
        plot(1:n, x_hist(:,3), 'b.-' );
        plot(1:n, agent.gt(:, 3), 'g.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);

        subplot(4, 1, 4); hold on; 
        plot(1:n, x_hist(:,4)/pi*180, 'b.-' );
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
            x_hist =  agents{veh_ids(i)}.my_x_hist();
            dxy = agents{veh_ids(i)}.gt(:, 1:2) - x_hist(:, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL MEAN: [%.3f] MEDIAN: [%.3f]\n', mean(dxy_all), median(dxy_all));
    
    pause(viz_frame_rate);
end

sol.agents = agents;
sol.system_setting = system_setting;
sol.problem = problem;
save(sprintf('solution_state_sharing_%i', simulation_scenario), 'sol')
%save(sprintf('solution_tracking_%i_%i', simulation_scenario, gps_error_tests(test_i)*10), 'solution_state_sharing')

end


