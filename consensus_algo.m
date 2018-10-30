
% remove nodes less than 3 connections
remove_idx = [];
for k1 = 1 : length(idx)
    if sum(A(k1, :)) < 3
        remove_idx = [remove_idx, k1];
    end
end

A(remove_idx, :) = [];
A(:, remove_idx) = [];
idx(remove_idx) = [];

% Prepare data
n_controls = -1;

coors = [];
for k1 = 1 : size(A, 1)
        xyz = epoch(idx(k1) == epoch(:, 2), 3:4);
        coors = [coors; xyz];
end
n_agents = size(coors, 1);

meas = [];
for k1 = 1 : size(A, 1)
    for k2 = 1 : size(A, 1)
        if (A(k1, k2) ~= 0)
            link = links(and( links(:, 2) == idx(k1), links(:, 3) == idx(k2)), 2:4);
            
            d = sqrt( (coors(k1, 1) - coors(k2, 1))^2 + (coors(k1, 2) - coors(k2, 2))^2 ) + normrnd(0, 0.20);
            meas = [meas; k1 k2 d];      
        end
    end
end

coors_init = coors;
idx_controls = 1:(n_agents-n_controls-1);
%coors_init(idx_controls, :) = coors_init(idx_controls, :) + (rand(length(idx_controls), size(coors, 2))-0.5)*1;
coors_init(idx_controls, :) = coors_init(idx_controls, :) + normrnd(0, 10, length(idx_controls), size(coors, 2));

x0l = coors_init';
x0l = x0l(:);
x0 = x0l; 

ci = (n_agents-n_controls):n_agents;     
x0( ((n_agents-n_controls)*2-1):n_agents*2 ) = [];

% Solve the problem
opts = optimoptions(@lsqnonlin,'Display','iter', 'Algorithm', 'levenberg-marquardt');
sol = lsqnonlin(@(x)obj_fn(x, x0l, meas, ci), x0, [], [], opts);

%opts = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
%sol = fminunc(@(x)norm(obj_fn(x, x0l, meas, ci)), x0, opts);

% Residuals
dr = obj_fn(sol, x0l, meas, ci);

%fprintf('Before: %.3f  After: %.3f\n', mean(obj_fn(x0, x0l, meas, ci)), mean(obj_fn(sol, x0l, meas, ci)));

coors_sol = reshape(sol, 2, length(x0)/2)';

figure(2); clf; hold on;
plot(coors(:, 1), coors(:, 2), 'g.', 'MarkerSize', 15);
plot(coors_init(:, 1), coors_init(:, 2), 'r.', 'MarkerSize', 15);
plot(coors_sol(:, 1), coors_sol(:, 2), 'b.', 'MarkerSize', 15);
plot(coors_sol(:, 1), coors_sol(:, 2), 'b+', 'MarkerSize', 10);
for k1 = 1 : size(meas, 1)
    plot([coors(meas(k1,1), 1), coors(meas(k1,2), 1)], [coors(meas(k1,1), 2), coors(meas(k1,2), 2)], 'r.-')
end
axis equal;
grid on;

% dx1 = norm( coors - coors_init );
% dx2 = norm( coors(1:size(coors_sol,1), :) - coors_sol );
dx1 = coors - coors_init;
dx1 = sqrt(sum((dx1(:)).^2) / length(dx1(:)));
%dx1 = mean(dx1(:));
dx2 = coors(1:size(coors_sol,1), :) - coors_sol;
dx2 = sqrt(sum((dx2(:)).^2) / length(dx2(:)));
%dx2 = mean(dx2(:));

title(sprintf("Before: %.3f After: %.3f", dx1, dx2));


return
%% Usage

% Some book keeping magic
x0l = coors';
x0l = x0l(:);
x0 = coors(1,:)'; % initial guess for unknown point

return;

Au = abs(A);
C = zeros(n_agents*2, n_agents*2);
for i = 1 : n_agents
    for j = (i+1) : n_agents
        if Au(i,j) == 1
            col = zeros(n_agents*2, 1);
            col((i-1)*2+1) = 1;
            col((i-1)*2+2) = 1;
            col((j-1)*2+1) = -1;
            col((j-1)*2+2) = -1;
            C = [C, col];
        end
    end
end

L = C*C';

n = 2*max(diag(L)); % max deg.
epsilon = 1/n/2;
epsilon = epsilon / 2;

%n = sum(diag(L)); 
%epsilon = 1/n;

%epsilon = 1e-9;
W = eye(size(L,1)) - epsilon*L; % weight matrix as graph's laplacian


