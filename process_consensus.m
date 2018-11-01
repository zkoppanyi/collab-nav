clear variables;

% Clear persistent variables
clear CoopAgent
clear CoopAgent2

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
viz_frame_rate = 0.1;                                                       % Frame rate for visualization
dt = 0.05;

% Error characterization
simulation_settings;
%system_setting.com_radius = 0.5;
sel_veh_id  = 17; 

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
x_m = zeros(n_veh*n_states, 1);
P_m = eye(n_veh*n_states, n_veh*n_states);

init_rewind = form_t;
init_i = 0;
init_time = 0.5/dt;
is_init = 0;
t = form_t;
iter_i = 0;

while 1

    if t >= to_t
        break;
    end
    
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
%     if isempty(find(veh_ids == sel_veh_id, 1))
%         continue;
%     end
    
    % Fixed number of vehicles
    rm_idx = [];
    for i = 1 : length(veh_ids)
        if isempty(find(veh_ids(i) == veh_ids_total))
            rm_idx = [rm_idx; i];
        end
    end
    veh_ids(rm_idx) = [];
    
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
        upt = [xy', v, b, t, is_init];  
        
        % Initialization
        if isempty(agents{veh_id})
           x_init = [xy; v; b];
           P_init = eye(4,4);  
           agents{veh_id} = CoopAgent(veh_id, x_init, P_init, t, n_veh, system_setting);
           %agents{veh_id} = CoopAgent2(veh_id, x_init, P_init, t, n_veh, system_setting);
           
            is_init == 1;
            init_i = 0; % apply intialzation
            init_rewind = t;
        end        
            
        % Build prediction: F, f, Q
        agents{veh_id}.build_predict();         
       
         % Measurement\observation        
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

            if ~isempty(agents{veh_id}.links)
                if isempty(find( agents{veh_id}.links(:, 2) == veh_id_j ))
                    agents{veh_id}.comm(agents{veh_id_j}, 'UWB'); 
                end
            else
                agents{veh_id}.comm(agents{veh_id_j}, 'UWB'); 
            end
            
            % Check that the other vehicle has already measured the this
            % vehicle
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
     
     % Initialize all vehciles: faster intialization
     % this step can be implemented as gossip algorithm
     if is_init
        for i = 1 : length(veh_ids) 
            for j = 1 : length(veh_ids) 
                agents{veh_ids(i)}.comm(agents{veh_ids(j)}, 'share-states');
            end
        end
     end
  
    %%
    % Consensus filter
    % Code based on:
    % https://ac.els-cdn.com/S0005109816300188/1-s2.0-S0005109816300188-main.pdf?_tid=d446a434-ef8e-48ca-a025-be3f6977a7a4&acdnat=1539389828_3e199ac4f616bfd57f1dbf569fdb3b61    
    
    [~, G, Adj] = get_local_conn(agents, veh_ids, veh_ids(1)); 
    %C = create_consensus_matrix(Adj, 'laplacian-const');
    C = create_consensus_matrix(Adj, 'max-degree');
    %C = create_consensus_matrix(Adj, 'laplacian-vary');
    
    for i = 1 : length(veh_ids)    

        agent = agents{veh_ids(i)};
        agent.build_predict(); 
        
        [h, H, z, R] = agent.build_update();                    
        z = z - afun_eval(h, agent.x) + H*agent.x;
        agent.dq     = H'*inv(R)*z;
        agent.dOmega = H'*inv(R)*H;
        
    end 
        
    % Simple average to test consensus
%     w = 1/length(veh_ids);
%     q_sum       = w * agents{veh_ids(1)}.q;
%     Omega_sum   = w * agents{veh_ids(1)}.Omega;
%     dq_sum      = w * agents{veh_ids(1)}.dq;
%     dOmega_sum  = w * agents{veh_ids(1)}.dOmega;
%     for i = 2 : length(veh_ids)
%         q_sum       = q_sum      + w * agents{veh_ids(i)}.q;
%         Omega_sum   = Omega_sum  + w * agents{veh_ids(i)}.Omega;
%         dq_sum      = dq_sum     + w * agents{veh_ids(i)}.dq;
%         dOmega_sum  = dOmega_sum + w * agents{veh_ids(i)}.dOmega;
%     end
%     for i = 1 : length(veh_ids)
%         agents{veh_ids(i)}.q        = q_sum;
%         agents{veh_ids(i)}.Omega    = Omega_sum;
%         agents{veh_ids(i)}.dq       = dq_sum;
%         agents{veh_ids(i)}.dOmega   = dOmega_sum;
%     end
    
    % Consensus
    Omegas_est = {};
    
    iter = [];
    %for ci = 1 : length(C)
    for ci = 1 : 50
        
         new_agents = agents;
         iter_col = zeros(length(veh_ids), 1);
         
         for i = 1 : length(veh_ids)
             
            agent = agents{veh_ids(i)};
            links = agent.links;  
            links = links(links(:, 2) < 1000, :);
            if isempty(links)
                continue
            end            
             
            % Save iterations
            idxs = agents{sel_veh_id}.idxs;
            iter_col(i)     = agent.q(idxs(1));
            iter_col(i)     = norm(agent.q);
            
            %w = C{ci}(i,i);
            w = C(i,i); 
            sum_w = w;
            q      = w * agent.q;
            Omega  = w * agent.Omega;            
            dq     = w * agent.dq;
            dOmega = w * agent.dOmega;            
                        
            for j = 1 : size(links, 1)
                idx = find(links(j, 2) == veh_ids);
                agent2 = agents{links(j, 2)};
                
                %w = C{ci}(i, idx);
                w      = C(i, idx); 
                sum_w = sum_w + w;
                q      = q       + w * agent2.q;
                Omega  = Omega   + w * agent2.Omega;
                dq     = dq      + w * agent2.dq;
                dOmega = dOmega  + w * agent2.dOmega;    
            end
                        
            agent.q         = q;
            agent.Omega     = Omega;
            agent.dq        = dq;
            agent.dOmega    = dOmega;         

            if sum(isnan(q)) ~= 0
                error('Somthing went wrong! NaN value in q!')
            end
            
            if abs(sum_w - 1) > 1e-5
                error('Somthing went wrong! sum_w ~= 1!')                
            end
            
            Omegas_est{ci} = inv(Omega);            
            agents{veh_ids(i)} = agent;        
            
         end         
         
         iter = [iter, iter_col];
         agents = new_agents;                     
    end
    
    figure(3); clf; hold on;
    plot(iter');
    title ('Convergence');
    
    
    %gamma = length(veh_ids);
    gamma = 1;
    Omegas = {};
    for i = 1 : length(veh_ids)    

        agent = agents{veh_ids(i)};
        
        agent.q     = agent.q       + gamma*agent.dq;
        agent.Omega = agent.Omega   + gamma*agent.dOmega;
        P = inv(agent.Omega);
        x = P*agent.q;
        
        agent.apply_update(x, P);
        agent.build_predict();  
        
        x_hat = afun_eval(agent.f, x);
        A = agent.F;
        W = inv(agent.Q);
        agent.Omega = W - W*A*inv(agent.Omega + A'*W*A)*A'*W;
        agent.q = agent.Omega*x_hat;
        Omegas{i} = inv(agent.Omega);        

        
        agents{veh_ids(i)} = agent;        
    end 

    %% Covariances
%     figure(3); clf; hold on; 
%     sidx = agents{sel_veh_id}.sidx;
%     for i = 1 : length(Omegas)
%         error_ellipse(Omegas{i}(sidx:(sidx+1), sidx:(sidx+1)), [0 0], 'style', 'b');
%         error_ellipse(Omegas_est{i}(sidx:(sidx+1), sidx:(sidx+1)), [0 0], 'style', 'r');
%     end
%     axis equal;
    

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
    
    figure(1); 
     if ~is_init
        title(sprintf("t= %.3f err=%.3f\n", t, mean(dxy_all)));
     end
    
    pause(viz_frame_rate);
end


sol.agents = agents;
sol.system_setting = system_setting;
sol.problem = problem;
save(['results\' sprintf('solution_consensus_%i', simulation_scenario)], 'sol')
%save(sprintf('solution_consensus_stat_%i_%i', simulation_scenario, system_setting.sigma_GPS*10), 'sol')

