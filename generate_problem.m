clear variables; 

%% Settings
const_bias_deg = 1.5;
const_bias_rad = const_bias_deg/180*pi;

max_connection = 10;
network_offset = [326819.60,4428871.90];
dt = 0.05;

%% Load SUMO files
folder = 'trajectories';
problem = {}; 
problem.vehicles = []; 
problem.roads = []; 
for k = 1 : 224
       
    fname = [folder '\vehicle_' num2str(k) '.txt'];
    s = dir(fname);
    fprintf('Vehicle #%i Size: %i bytes\n', k, s.bytes);
    
    if s.bytes == 0
        continue;
    end
    
    data = dlmread(fname);   
    %data(:,2:3) = data(:,2:3) + network_offset;
                
    db = [0; diff(data(:,4))];
    dt = data(2, 1) - data(1, 1);

%         % Make it Ackerman steering vehicle
%         v = problem.vehicles{k}(:,5);
%         delta_2 = zeros(length(v), 1);
%         idx = find(v ~= 0);
%         delta_2(idx) = atan( (L * 1./v(idx)) .* (db(idx)/180*pi/dt) );
%         problem.vehicles{k}(:,4) = problem.vehicles{k}(:,4) + delta_2/pi*180;
%         %db = [0; diff(problem.vehicles{k}(:,4))];

    %bias = (problem.vehicles{k}(:,1) - problem.vehicles{k}(1,1))*0.01;
    data = [data, db/dt + const_bias_deg];       

    t1 = data(:,1);
    t2 = min(t1):dt:max(t1);
    data2 = zeros(length(t2), size(data, 2));
    for i = 1 : size(data, 2)
        data2(:, i) = interp1(t1, data(:,i), t2, 'spline');
    end
   
    problem.vehicles = [problem.vehicles; data2(:, 1), repmat(k, size(data2,1), 1), data2(:, 2:end)]; 



end

% Get the roads for visualization
problem.roads = problem.vehicles(:, 3:4);

save('problem', 'problem');

return;
%% Further processing: creating rel_obs and epoch structures
% Data structures
% -------------------
% 
% 1., problem.epoch{idx}
% time [s], vehicle id, X, Y, angle, speed, (list: (connected vehicle, distance)
%
% 2., problem.rel_obs
% source vehicle id, target vehicle id, distance
%
% 3., problem.vehicles{idx}
% time, X, Y, angle, speed
%

load('problem3')



% Get observations
problem.rel_obs = [];
problem.epochs = {};
iter = 0;
llh_all = [];
for ti = 0 : 0.2 : 1000
    
    iter = iter + 1;
    problem.epochs{iter} = [];    
    figure(1); clf; hold on;
   
    plot(roads(:, 1), roads(:, 2), 'k.');
    n_veh = 0;
    xyz = [];
    
    for k = 1 : length(problem.vehicles)
    %for k = 1 : 10
        if size(problem.vehicles{k}, 1) ~= 0
            idx = find( abs(problem.vehicles{k}(:, 1) - ti) < 0.0001 );   
            if ~isempty(idx)
                
                new_line = zeros(length(idx), size(problem.vehicles{k}, 2)+1);
                new_line(:, 1:size(problem.vehicles{k}, 2) - 1) = problem.vehicles{k}(idx, 2:end);
                utm_xyz = new_line(:,1:2) + network_offset;
                [lat, lon] = utm2deg(utm_xyz(1),utm_xyz(2),'17 N'); llh_all = [llh_all; lat lon];
                new_line(:, size(problem.vehicles{k}, 2):end) = [lat lon];                
                xyz = [xyz; k, new_line];
                
                plot(problem.vehicles{k}(idx, 2), problem.vehicles{k}(idx, 3), 'r.', 'MarkerSize', 15);
                n_veh = n_veh + 1;                                
            end
        end
    end
    
    %fprintf('%.1f %i\n', ti, size(traj, 1));
    
    if ~isempty(xyz)
        
        %ids = unique(traj(:, 1));
        for i = 1 : size(xyz, 1)
            d = sqrt((xyz(i, 2) - xyz(:, 2)).^2 + (xyz(i, 3) - xyz(:, 3)).^2);
            idx = find(and(d < 200, d ~= 0));

            obs = zeros(1, max_connection*2);
            for j = 1 : min(max_connection, length(idx))
                plot( [ xyz(i, 2), xyz(idx(j), 2) ], [ xyz(i, 3), xyz(idx(j), 3) ] , 'b.-' )
                
                d_meas = d(idx(j)); % + normrnd(0, 0.5);
                problem.rel_obs = [problem.rel_obs; ti xyz(i, 1) xyz(idx(j), 1) d_meas];
                
                obs(2*j-1) =  xyz(idx(j), 1);
                obs(2*j) =  d_meas;                
            end
            problem.epochs{iter} = [problem.epochs{iter}; ti xyz(i, :) obs];

        end
    
    end
       
    
    title(sprintf('time: %.1f; # of veh.: %i', ti, n_veh));
    axis equal;
    grid on;
    
    pause(0.01);
end



