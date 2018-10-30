clear variables;

% Papers about vehicle models
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5472738
% https://benthamopen.com/contents/pdf/TOASJ/TOASJ-3-13.pdf

load('problem3');

%% Settings 

dt          = 0.2;                   % Sampling rate
AOI         = [750 1100; 1250 1550]; % Area of interest
sel_veh_id  = 13;                    % Selected vehicle ID
viz_frame_rate = 0.3;                % Frame rate for visualization
L = 1;

% Error characteristics
sigma_GPS = 1.0;
sigma_v = 0.1;
sigma_IMU = (0.1/180*pi);
const_bias_deg = 1.5;
const_bias_rad = const_bias_deg/180*pi;

%% Initialization
roads = problem.roads;

for k = 1 : length(problem.vehicles)
    
    % Add bias to heading
    if ~isempty( problem.vehicles{k})
        
        db = [0; diff(problem.vehicles{k}(:,4))];
        
        % Make it Ackerman steering vehicle
        v = problem.vehicles{k}(:,5);
        delta_2 = zeros(length(v), 1);
        idx = find(v ~= 0);
        delta_2(idx) = atan( (L * 1./v(idx)) .* (db(idx)/180*pi/dt) );
        problem.vehicles{k}(:,4) = problem.vehicles{k}(:,4) + delta_2/pi*180;
        %db = [0; diff(problem.vehicles{k}(:,4))];
        
        %bias = (problem.vehicles{k}(:,1) - problem.vehicles{k}(1,1))*0.01;
        problem.vehicles{k} = [problem.vehicles{k}, db/dt + const_bias_deg];       
        
    end
    
end

% Get timestamps
ts = problem.rel_obs(:,1);
ts = unique(ts);

figure(1);
set(gcf,'Position',[50 200 500 500])
figure(2); clf; hold on;                
set(gcf,'Position',[600 200 800 700])

%% Algorithm
agents{length(problem.vehicles)} = {};
dxy_all = [];

%for k = 250 : length(problem.epochs)
for k = 500 : 800

    epoch = problem.epochs{k};
    if isempty(epoch)
        continue;
    end
    % Find trajectory within the AOI
    epoch = epoch(and(and(and(epoch(:,3) > AOI(1,1), epoch(:,3) < AOI(1,2)), epoch(:,4) > AOI(2,1)), epoch(:,4) < AOI(2,2)), :);
    t = epoch(1,1);
    
    %fprintf("t= %.1f\n", t);   
    
    % Vehicles at the current epoch
    veh_ids = unique(epoch(:,2));         
    
    if isempty(find(veh_ids == sel_veh_id, 1))
        continue;
    end
    
    figure(1); clf; hold on;
    plot(roads(:, 1), roads(:, 2), 'k.');

    for i = 1 : length(veh_ids)
        
        veh_id = veh_ids(i);
        
        %if veh_id ~= sel_veh_id, continue; end
        %if floor(veh_id) ~= veh_id, continue; end
        
        % Get internal observations
        veh = problem.vehicles{veh_id};
        [val, idx] = min(abs(veh(:,1) - t));
        if val < 0.01
            inter_obs = veh(idx,:);
        else
            %inter_obs = [];
            continue;
        end
                 
        % Get internal observation
        xy = inter_obs(1,2:3)';
        b  = inter_obs(4)/180*pi;
        v  = inter_obs(5);
        db = inter_obs(6)/180*pi;
        gt = [xy', v, b, db - const_bias_rad, const_bias_rad, 0.5, 0.5];              
        
        % Initialization
        if isempty(agents{veh_id})
           agents{veh_id}.x = [xy; v; b; db; const_bias_rad+1/180*pi; 0.5; 0.5];
           P = eye(8,8);  P(6,6) = 1; 
           agents{veh_id}.P = P;
           agents{veh_id}.x_all = [];
           agents{veh_id}.z = [];
           agents{veh_id}.gt_all = []; 
           agents{veh_id}.t0 = t; 
           agents{veh_id}.z_gps = []; 
        end                

        % Collect data and form matrices
        x = agents{veh_id}.x;
        P = agents{veh_id}.P;      
        F = eye(8, 8);
        F(1,3) = sin(x(4))*dt;
        F(2,3) = cos(x(4))*dt;
        F(4,5) = dt;
        %F(7,7) = -1/dt;
        %F(8,8) = -1/dt;
        
        %F = [1 0 sin(x(4))*dt 0 0 0; 0 1 cos(x(4))*dt 0 0 0; 0 0 1 0 0 0; 0 0 0 1 dt 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
        Q = diag([1, 1, 0.1, (0.1/180*pi)^2, (0.1/180*pi)^2, (0.0001/180*pi)^2 1 1]);
        R = diag([sigma_GPS^2, sigma_GPS^2, sigma_v^2, sigma_IMU^2, sigma_IMU^2, sigma_IMU^2, 1 1]);
        
        % State predication
        x = F*x;
        P = F*P*F' + Q;
        
        % Updates
        if mod(k, 5) == 0
            H = zeros(8,8); H(1,1) = 1; H(2,2) = 1; 
            xy_meas = xy + normrnd(0, sigma_GPS, 2, 1);
            z = [xy_meas; 0; 0; 0; 0; 0; 0];
            
            % Bearing update from the GPS
            if ~isempty(agents{veh_id}.z_gps)
                z_prev = agents{veh_id}.z_gps(end,:);
                b_gps = sign(b)*abs(atan2(z(1)-z_prev(1), z(2)-z_prev(2)));
                if abs(b_gps - x(4)) < 10/180*pi
                    H(4, 4) = 1;
                    z(4) = b_gps;
                    R(4,4) = 10;
                    %fprintf('Bearing update: %.1f, %.1f\n', b, b_gps);    
                end
            end
            agents{veh_id}.z_gps = [agents{veh_id}.z_gps; xy_meas']; 
        else
            H = zeros(8,8); H(3,3) = 1; H(5,5) = 1; H(5,6) = 1; 
            z = [0; 0; v + normrnd(0, sigma_v); 0; db + normrnd(0, sigma_IMU)/180*pi; 0; 0; 0]; 
        end
        
        y = z - H*x;
        S = R + H*P*H';
        K = P*H'*inv(S);
        x = x + K*y;
        P = (eye(size(K, 1), size(K, 1)) - K*H)*P*(eye(size(K, 1), size(K, 1)) - K*H)'+K*R*K';
        y = z - H*x;
        
        % Filter reset, if (1) estimated bias is far from ground truth and
        % (2) the filter is alive more than 2 cycles
        if and(abs(x(6) - const_bias_rad) > 1/180*pi, abs(agents{veh_id}.t0 - t) > dt*5)
            agents{veh_id} = {};
            fprintf("Reset: %i\n", veh_id);
            continue;
        end
               
        % Save results
        agents{veh_id}.x = x;
        agents{veh_id}.P = P;        
        agents{veh_id}.x_all = [agents{veh_id}.x_all; x'];
        agents{veh_id}.gt_all = [agents{veh_id}.gt_all; gt];     
        agents{veh_id}.z = [agents{veh_id}.z; z'];        

        % Visualization
        plot(x(1), x(2), 'b.', 'MarkerSize', 15);
        plot(xy(1), xy(2), 'g.', 'MarkerSize', 15);  
        plot(z(1), z(2), 'r.', 'MarkerSize', 15); 
        
        if veh_id == sel_veh_id
            text(x(1), x(2)+20, sprintf('%i', veh_id), 'BackgroundColor', 'w', 'Color', 'r');
        else
            text(x(1), x(2)+20, sprintf('%i', veh_id), 'BackgroundColor', 'w');
        end
        
%         if veh_id == sel_veh_id
%             % Observaibility
%             % See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.30.2399&rep=rep1&type=pdf
%             O = H; 
%             for l = 1 : 10
%                 O = [O; H*F^l];
%             end
%             fprintf('O = %i \n', rank(O));
%         end

    end    
        
    grid on; axis equal;    
    xlim(AOI(1,1:2)); ylim(AOI(2,1:2));       
    set(gca, 'FontSize', 12);    
    xlabel('[m]'); ylabel('[m]');
    
    if ~isempty(agents{sel_veh_id})
        
        n = size(agents{sel_veh_id}.x_all, 1);
        dxy = (agents{sel_veh_id}.x_all(end, :)-agents{sel_veh_id}.gt_all(end, :))';

        figure(2); clf; hold on;                
        %set(gcf,'Position',[600 200 800 700])
         
        subplot(6, 1, 1); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 1), 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 1), 'g.-' );
        title(sprintf('[%i] X= %.3f', sel_veh_id, dxy(1)));
        set(gca, 'FontSize', 10);

        subplot(6, 1, 2); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 2), 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 2), 'g.-' );
        title(sprintf('Y= %.3f', dxy(2)));
        set(gca, 'FontSize', 10);

        subplot(6, 1, 3); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 3), 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 3), 'g.-' );
        title(sprintf('v= %.3f', dxy(3)));
        set(gca, 'FontSize', 10);
        
        subplot(6, 1, 4); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 4)/pi*180, 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 4)/pi*180, 'g.-' );
        title(sprintf('bearing= %.3f', dxy(4)/pi*180));  
        set(gca, 'FontSize', 10);
        
        subplot(6, 1, 5); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 5)/pi*180, 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 5)/pi*180, 'g.-' );
        title(sprintf('\\Delta bearing= %.3f', dxy(5)/pi*180));
        set(gca, 'FontSize', 10);
        
        subplot(6, 1, 6); hold on; 
        plot(1:n, agents{sel_veh_id}.x_all(:, 6)/pi*180, 'b.-' );
        plot(1:n, agents{sel_veh_id}.gt_all(:, 6)/pi*180, 'g.-' );
        title(sprintf('\\Delta bearing bias = %.3f', dxy(6)/pi*180)); 
        set(gca, 'FontSize', 10);
        
    end
    
    %clc   
    fprintf('Time: %.1f\n', t);    
    for i = 1 : length(veh_ids)
        if ~isempty(agents{veh_ids(i)})
            dxy = agents{veh_ids(i)}.gt_all(:, 1:2) - agents{veh_ids(i)}.x_all(:, 1:2);
            dxy = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
            xy_rmse = sqrt(sum(dxy.^2) / length(dxy));
            xy_std = std(dxy);
            fprintf('#%i RMSE: [%.3f]  STD: [%.3f]\n', veh_ids(i), xy_rmse, xy_std);
            dxy_all  = [dxy_all; dxy(end)];
        end
    end
    fprintf('TOTAL RMSE: [%.3f]  MEAN: [%.3f] MEDIAN: [%.3f]\n', sqrt(sum(dxy_all.^2) / length(dxy_all)), mean(dxy_all), median(dxy_all));
    
    
    pause(viz_frame_rate);

end



