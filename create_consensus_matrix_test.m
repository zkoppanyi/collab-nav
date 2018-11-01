clear variables; clc;

n_agents = 150;
xv_gt = normrnd(15,150,n_agents,1);
iters = xv_gt;
xv = xv_gt;

% Generate connected incidence matrix
n_max_conn = 5;
n_min_conn = 3;
A = zeros(n_agents, n_agents);
for i = 1 : n_agents 
    conn = randperm(n_agents);
    conn = conn(conn ~= i);
    n_conn = sum(A(i, :));
    j = 0;
    while and(n_conn < n_min_conn, j < length(conn))
        j = j + 1;
        if and(A(i, conn(j)) == 0, sum(A(conn(j), :)) < n_max_conn)
            A(i, conn(j)) = 1;
            A(conn(j), i) = 1;
            n_conn = n_conn + 1;
        end
    end    
end

figure(1);
G = graph(A);
plot(G);

%W = create_consensus_matrix(A, 'laplacian-const');    
%W = create_consensus_matrix(A, 'max-degree');    
W = create_consensus_matrix(A, 'metropolis');    
for k = 1 : 200
    xv = W*xv;
    iters = [iters, xv];
end

% W = create_consensus_matrix(A, 'laplacian-vary');
% for k = 1 : length(W)
%     xv = W{k}*xv;
%     iters = [iters, xv];
% end

fprintf('%.3f =? %.3f [std= %.6f max= %.6f] \n', mean(xv_gt), iters(1, end), std(iters(:, end)), max(abs(iters(:, end) - mean(xv_gt))) );

%iterations
% Au = abs(A);
% for k = 1 : 60    
%     xv2 = xv;
%     for i = 1 : n_agents
%         sxv = W(i,i)*xv2(i);
%         for j = 1 : n_agents
%             if and(Au(i,j) == 1, i~=j)
%                 sxv = sxv + W(i,j)*xv2(j);
%             end;
%         end;
%         xv(i) = sxv ;
%     end;
%     iters = [iters, xv];
% end;
% xv

figure(2);
plot(iters');