function W = create_consensus_matrix(A, type, is_check)

    if nargin < 3
        is_check = 1;        
    end
    
    if size(A, 1) ~= size(A, 2)
        error('Adjacency matrix is not square!');
    end
    
    if sum(sum(A - A')) ~= 0
        error('Adjacency matrix is not symmetric!');
    end

    my_eps                  = 1e-5;
    needed_right_stoch      = 0;
    needed_left_stoch       = 0;
    needed_symmetry         = 0;
    
    if is_check
                  
         % Using Matlab's function
         G = graph(A);
         v = conncomp(G);
         if max(v) ~= 1
            warning('Graph is not connected!');
         end
        
         L = compute_laplacian(A);
         [~, S, ~] = eig(L);
         v = diag(S);
         if v(2) < my_eps
            warning('Graph is not connected!');
        end
    end
    
    
    if strcmp(type, 'laplacian-const')
        
        % Constant edge weights
        
        L = compute_laplacian(A);
        [~, S, ~] = eig(L);
        nzv = diag(S);
        nzv=nzv(nzv>my_eps); % non zero eigenvalues
        nzv = sort(nzv, 'descend'); % this is not needed just for nicer convergence figures
        alpha = 2 / (min(nzv) + max(nzv)); 
        W = eye(size(L,1), size(L,2)) - alpha*L;
        
        needed_right_stoch      = 1;
        needed_left_stoch       = 1;
        needed_symmetry         = 1;        
        
    elseif strcmp(type, 'laplacian-vary')    
        
        % Based on https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6315398
        % Allows reaching the exact average consensus in a finite number of 
        % steps when perfect data exchange is considered
        
        L = compute_laplacian(A);
        [~, S, ~] = eig(L);
        nzv = diag(S);
        nzv = nzv(nzv > my_eps); % non zero eigenvalues
        nzv = sort(nzv, 'descend'); % this is not needed just for nicer convergence figures        
        
        n_agents = size(A, 1);
        xv_gt = normrnd(15,150,n_agents,1);
        xv = xv_gt;
        
        W = {};
        for k = 1 : length(nzv)
             W_temp = eye(size(L,1), size(L,2)) - 1/nzv(k)*L;   
             
             % checking for round-off error
             xv = W_temp*xv;
             dxv = norm(xv - xv_gt);
             if  dxv > 1/my_eps
                 warning(sprintf('Poor convergence, possibly round-off error\n Exiting... Weights might not be correct! detected difference: %.6f', dxv));
                 break;
             end
        
             W{k} = W_temp;
        end
        
        needed_right_stoch      = 1;
        needed_left_stoch       = 1;    
        needed_symmetry         = 1;
        
    elseif strcmp(type, 'max-degree')    
        
        % See: http://web.stanford.edu/~boyd/papers/pdf/lms_consensus.pdf
        % Uses only local information, does not require the network
        % topology
        
        d = -Inf;
        for i = 1 : size(A, 1)
            di = sum(A(i, :));
            if di > d
                d = di;
            end
        end
        
        n_agents = size(A, 1);
        W = zeros(n_agents, n_agents);
        for i = 1 : size(A, 1)
            for j = 1 : size(A, 1)  
                if i == j
                    di = sum(A(i, :));
                    W(i,j) = 1 - di / (d + 1);
                elseif A(i,j) == 1                    
                    W(i,j) = 1 / (d + 1);
                end
            end
        end
        
        needed_right_stoch      = 1;
        needed_left_stoch       = 1;  
        needed_symmetry         = 1;
        
   elseif strcmp(type, 'metropolis')    
        
        % See: http://web.stanford.edu/~boyd/papers/pdf/lms_consensus.pdf
        % Uses only local information, does not require the network
        % topology
        % Does require the graph to be bipirate
        
        d = -Inf;
        for i = 1 : size(A, 1)
            di = sum(A(i, :));
            if di > d
                d = di;
            end
        end
        
        n_agents = size(A, 1);
        W = zeros(n_agents, n_agents);
        
        for i = 1 : n_agents
            di = sum(A(i,:));    
            for j = 1 : n_agents
               if i ~= j
                   if A(i,j) == 1
                       dj = sum(A(j,:));  
                       W(i,j) = min(1/di, 1/dj);
                   else
                       W(i,j) = 0;
                   end
               else
                   for k = 1 : n_agents
                       if A(i,k) == 1
                           dk = sum(A(k,:));  
                           W(i,i) = W(i,i) + max(0, 1/di-1/dk);
                       end
                   end
               end

            end
        end

        needed_right_stoch      = 1;
        needed_left_stoch       = 1;  
        needed_symmetry         = 1;        
        
    end
    
    % Checking... 
    if is_check
        
        Wa = {};
        if ~iscell(W)
            Wa{1} = W;
        else
            Wa = W;
        end

        for k = 1 : length(Wa)

            % Is it left stochastic?
            r = Wa{k}*ones(size(Wa{k}, 1), 1);
            if and( abs(sum(r) - size(Wa{k}, 1)) > my_eps, needed_left_stoch)
                error('Result is not left stochastic!');
            end

            % Is it right stochastic?
            r = ones(1, size(Wa{k}, 1)) * Wa{k};
            if and( abs(sum(r) - size(Wa{k}, 1)) > my_eps, needed_right_stoch)
                error('Result is not right stochastic!');
            end

            if and(sum(sum(Wa{k} - Wa{k}')) ~= 0, needed_symmetry)
                error('Result is not symmetric!');
            end
           
        end      
    end
    
end
    
function L = compute_laplacian(A)
    n_agents = size(A, 1);
    C = [];
    for i = 1 : n_agents
        for j = (i+1) : n_agents
            if A(i,j) ~= 0
                col = zeros(n_agents, 1);
                col(i) = 1;
                col(j) = -1;
                C = [C, col];
            end
        end
    end
    L = C*C';
end
