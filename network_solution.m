clear variables;

load('problem2');

roads = [];
for k = 1 : length(problem.vehicles)
    if size(problem.vehicles{k}, 1) ~= 0
        roads = [roads; problem.vehicles{k}(:, 2), problem.vehicles{k}(:, 3)];
    end
end

ts = problem.rel_obs(:,1);
ts = unique(ts);

% [pX pY] = meshgrid( min(roads(:,1)):50:max(roads(:,1)), min(roads(:,2)):50:max(roads(:,2)));
% pX = pX(:); pY = pY(:); 
% static_net = [];
% for k = 1 : length(pX)
%     [val, idx] = min( (pX(k) - roads(:,1)).^2 + (pY(k) - roads(:,2)).^2 );
%     if val < 50
%         static_net  = [static_net; pX(k) pY(k) ];
%     end
% end

%for k = 250 : length(problem.epochs)
for k = 500 : 500

    epoch = problem.epochs{k};
    if isempty(epoch)
        continue;
    end
    t = epoch(1,1);
    
    links = problem.rel_obs(  problem.rel_obs(:, 1) == t, :);
    n = size(links, 1);
    fprintf("t= %.1f n= %i\n", t, size(links, 1) );   
    
    
    if n == 0
        continue
    end
    
    figure(1); clf; hold on;
    plot(roads(:, 1), roads(:, 2), 'k.');
    plot(epoch(:, 3), epoch(:, 4), 'b.', 'MarkerSize', 15);
    %plot(static_net(:, 1), static_net(:, 2), 'r.', 'MarkerSize', 15);
    grid on;
    axis equal;
    
    for i = 1 : size(links, 1)
       idx1 = find(links(i, 2) == epoch(:, 2));
       idx2 = find(links(i, 3) == epoch(:, 2));                  
       
       plot([epoch(idx1, 3) epoch(idx2, 3)], [epoch(idx1, 4) epoch(idx2, 4)], 'b-', 'LineWidth', 3);       
       
    end
        
    G = graph(links(:, 2), links(:, 3));    
    [bin, binsize] = conncomp(G);
    bidx = find(binsize >= 3);
    fprintf("\t n: %i (", length(bidx));   
    
    net = {};
    nidx = {};
    for j = 1 : length(bidx)
        nidx{j} = find(bin == bidx(j));
        sub_network = subgraph(G, nidx{j});
        net{j} = full(adjacency(sub_network));
        fprintf("%i, ", size(net, 1) );   
        %figure(j); clf; hold on; plot(SG);
    end
    fprintf(")\n"); 
   
    if isempty(net)
        continue
    end
    
    A = net{1};
    idx = nidx{1};
    
    consensus_algo
        
    pause(1);
end



