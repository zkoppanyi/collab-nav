classdef CoopAgent < matlab.mixin.Copyable
    properties
        
        id              = -1;
        sid             = -1;
        max_neighbors   = -1;

        n_states = 0;
        x        = [];
        P        = [];        

        z_int    = [];
        z_ext    = [];

        x_hist    = [];        
        z_int_hist   = [];
        z_ext_hist   = [];

        gt       = [];                                                      % ground truth
        t0       = 0;  

        F        = [];
        f        = {};
        H_int    = [];
        h_int    = {};
        R_int    = [];
        H_ext    = [];
        h_ext    = {};
        Q        = [];
        K        = [];
        
        q        = [];
        Omega    = [];
        dq        = [];
        dOmega    = [];
        
        links    = [];
        tracked  = [];                                                      % tracked agents
        
        system_setting = {};
        
        sidx = -1;
        idxs = [];
        
    end
            
    methods        
        
        function obj = CoopAgent(id, x_init, P_init, t_init, max_neighbors, system_setting)
            
            if id ~= -1
                persistent nsidx
                if isempty(nsidx)
                   nsidx = 1;
                else
                    nsidx = nsidx + 1;
                end
                obj.sid             = nsidx;
            end
            
            obj.id              = id;            
            obj.max_neighbors   = max_neighbors;
            %obj.n_states        = length(x_init);
            obj.n_states        = 4;
            obj.x               = zeros(max_neighbors*obj.n_states, 1);
            obj.P               = eye(length(obj.x));
            obj.t0              = t_init;
            obj.system_setting  = system_setting;
            
            obj.q               = zeros(max_neighbors*obj.n_states, 1);
            obj.Omega           = eye(length(obj.q), length(obj.q));
            
            if id ~= -1
                obj.sidx = (obj.sid-1)*obj.n_states+1;
                obj.idxs = obj.sidx:(obj.sidx+obj.n_states-1);   
                obj.x(obj.idxs) = x_init;
                obj.q(obj.idxs) = x_init;
            end
            
        end
        
        function st = my_x_hist(obj)
            st = obj.x_hist(:, obj.idxs);
        end
        
        function apply_update(obj, x, P)
            
                if (size(x, 1) > size(x, 2))
                    x2 = x';
                end
            
                obj.x = x;
                obj.P = P;
                obj.Omega = inv(P);
                obj.q = obj.Omega * x;
                
                obj.x_hist = [obj.x_hist; x2, P(:)'];
        end
            
        function build_predict(obj)
           
            n_total = obj.max_neighbors * obj.n_states;            
            obj.F = zeros(n_total, n_total);
            obj.Q = zeros(n_total, n_total);
            obj.f = cell(n_total, 1);
            %obj.f = afun_create(n_total, @(x) 0)
            
            for j = 1 : obj.max_neighbors
                
                lsidx = (j - 1) * obj.n_states + 1;
                lidxs = lsidx:(lsidx + obj.n_states - 1); 
                obj.F(lidxs,lidxs) = [1 0 sin(obj.x(lidxs(4)))*obj.system_setting.dt 0; 0 1 cos(obj.x(lidxs(4)))*obj.system_setting.dt 0; 0 0 1 0; 0 0 0 1];

                obj.f{lsidx+0} = @(x) x(lidxs(1)) +  x(lidxs(3))*sin(x(lidxs(4))) * obj.system_setting.dt;
                obj.f{lsidx+1} = @(x) x(lidxs(2)) +  x(lidxs(3))*cos(x(lidxs(4))) * obj.system_setting.dt;
                obj.f{lsidx+2} = @(x) x(lidxs(3));
                obj.f{lsidx+3} = @(x) x(lidxs(4));

                obj.Q(lidxs,lidxs) = diag([0.5, 0.5, 0.05, (2.5/180*pi)^2]);
                
            end
        end
      
        function build_int_update(obj, type, upt)
            
            if (size(upt, 1) > size(upt, 2))
                upt = upt';
            end
            
            obj.gt = [obj.gt; upt];
        
            lsidx = obj.sidx;
            
            if or(or(or( strcmp(type, 'GPS'), strcmp(type, 'IMU')), strcmp(type, 'GPS+IMU')), strcmp(type, 'ZUPT'))
                n_total = obj.n_states*obj.max_neighbors;
                obj.H_int = zeros(n_total, n_total);
                obj.R_int = eye(n_total, n_total);
                obj.z_int = zeros(n_total,1);
                obj.h_int = afun_create(n_total, @(x) 0);
            end

            if strcmp(type, 'GPS')

                obj.H_int(lsidx+0,lsidx+0)  = 1; 
                obj.H_int(lsidx+1,lsidx+1)  = 1; 
                
                obj.R_int(lsidx,lsidx)     = obj.system_setting.sigma_GPS^2; 
                obj.R_int(lsidx+1,lsidx+1) = obj.system_setting.sigma_GPS^2;

                obj.h_int{lsidx+0}   =  @(x) x(lsidx);
                obj.h_int{lsidx+1}   =  @(x) x(lsidx+1);        
                obj.h_int{lsidx+2}   =  @(x) 0;
                obj.h_int{lsidx+3}   =  @(x) 0;

                obj.z_int(lsidx+0)     =  upt(1) + normrnd(0, obj.system_setting.sigma_GPS);        
                obj.z_int(lsidx+1)     =  upt(2) + normrnd(0, obj.system_setting.sigma_GPS);                 
                obj.z_int(lsidx+2)     =  0;        
                obj.z_int(lsidx+3)     =  0; 

            elseif strcmp(type, 'IMU')

                obj.H_int(lsidx+2,lsidx+2) = 1; 
                obj.H_int(lsidx+3,lsidx+3) = 1; 
                
                obj.R_int(lsidx+2,lsidx+2) = obj.system_setting.sigma_v^2; 
                obj.R_int(lsidx+3,lsidx+3) = obj.system_setting.sigma_IMU^2;

                obj.h_int{lsidx+0}     =  @(x) 0;
                obj.h_int{lsidx+1}     =  @(x) 0;
                obj.h_int{lsidx+2}     =  @(x) x(lsidx+2);
                obj.h_int{lsidx+3}     =  @(x) x(lsidx+3);

                obj.z_int(lsidx+0)     =  0;        
                obj.z_int(lsidx+1)     =  0; 
                obj.z_int(lsidx+2)     =  upt(3) + normrnd(0, obj.system_setting.sigma_v); 
                obj.z_int(lsidx+3)     =  upt(4) + normrnd(0, obj.system_setting.sigma_IMU);  

            elseif strcmp(type, 'GPS+IMU')

                obj.H_int(lsidx+0,lsidx+0)  = 1; 
                obj.H_int(lsidx+1,lsidx+1)  = 1; 
                obj.H_int(lsidx+2,lsidx+2)  = 1; 
                obj.H_int(lsidx+3,lsidx+3)  = 1;
                
                obj.R_int(lsidx+0, lsidx+0) = obj.system_setting.sigma_GPS^2; 
                obj.R_int(lsidx+1, lsidx+1) = obj.system_setting.sigma_GPS^2;
                obj.R_int(lsidx+2, lsidx+2) = obj.system_setting.sigma_v^2; 
                obj.R_int(lsidx+3, lsidx+3) = obj.system_setting.sigma_IMU^2;

                obj.h_int{lsidx+0}     =  @(x) x(lsidx+0);
                obj.h_int{lsidx+1}     =  @(x) x(lsidx+1);
                obj.h_int{lsidx+2}     =  @(x) x(lsidx+2);
                obj.h_int{lsidx+3}     =  @(x) x(lsidx+3);

                obj.z_int(lsidx+0)     =  upt(1) + normrnd(0, obj.system_setting.sigma_GPS);        
                obj.z_int(lsidx+1)     =  upt(2) + normrnd(0, obj.system_setting.sigma_GPS);    
                obj.z_int(lsidx+2)     =  upt(3) + normrnd(0, obj.system_setting.sigma_v); 
                obj.z_int(lsidx+3)     =  upt(4) + normrnd(0, obj.system_setting.sigma_IMU);      
                
            elseif strcmp(type, 'ZUPT')
                    
                obj.H_int(lsidx+0,lsidx+0)  = 1; 
                obj.H_int(lsidx+1,lsidx+1)  = 1; 
                obj.H_int(lsidx+2,lsidx+2)  = 1; 
                obj.H_int(lsidx+3,lsidx+3)  = 1;
                
                obj.R_int(lsidx+0, lsidx+0) = 0.1^2; 
                obj.R_int(lsidx+1, lsidx+1) = 0.1^2;
                obj.R_int(lsidx+2, lsidx+2) = 0.01^2; 
                obj.R_int(lsidx+3, lsidx+3) = 0.01^2;

                obj.h_int{lsidx+0}     =  @(x) x(lsidx+0);
                obj.h_int{lsidx+1}     =  @(x) x(lsidx+1);
                obj.h_int{lsidx+2}     =  @(x) x(lsidx+2);
                obj.h_int{lsidx+3}     =  @(x) x(lsidx+3);

                obj.z_int(lsidx+0)     =  upt(1);        
                obj.z_int(lsidx+1)     =  upt(2);    
                obj.z_int(lsidx+2)     =  upt(3);  
                obj.z_int(lsidx+3)     =  upt(4);          

            end
            
            % Save some history data             
            obj.z_int_hist = [obj.z_int_hist; obj.z_int'];     
            
        end
        
        function [h, H, z, R] = build_update(obj)
            H = [obj.H_int; obj.H_ext]; 
            z = [obj.z_int; obj.z_ext]; 
            R_add = diag( ones(length(obj.z_ext), 1)*obj.system_setting.sigma_UWB^2 );
            R = obj.R_int;
            R = [R, zeros(size(R, 1), size(R_add, 2));  zeros(size(R_add, 1), size(R, 1)), R_add]; 
            h = afun_concat(obj.h_int, obj.h_ext);
        end
        
        function comm(obj, agent2, type, varargin)           
          
          if strcmp(type, 'UWB') 
            
            idxs_j = agent2.idxs;
            obj.x(idxs_j) = agent2.x(idxs_j);
            obj.P(idxs_j, idxs_j) = agent2.P(idxs_j, idxs_j);
            obj.tracked = [obj.tracked; agent2.id];
          
            xy_gt_i   = obj.gt(end, :);
            xy_gt_j   = agent2.gt(end, :);
            
            % This means that the vehicle left AOI, so the coordinates are
            % not updated; skip this communication link
            if ~(and(and(and(xy_gt_j(1) > obj.system_setting.AOI(1,1), ...
                    xy_gt_j(1) < obj.system_setting.AOI(1,2)), ...
                    xy_gt_j(2) > obj.system_setting.AOI(2,1)), ...
                    xy_gt_j(2) < obj.system_setting.AOI(2,2)))
                
                return;
            end
            
            obj.links = [obj.links; obj.id, agent2.id]; 
            d = sqrt((xy_gt_i(1) - xy_gt_j(1))^2 + (xy_gt_i(2) - xy_gt_j(2))^2);
            d = d + normrnd(0, obj.system_setting.sigma_UWB);                  % Simulate ranging      
            
            if d > obj.system_setting.com_radius
               return;
            end
            
            sidx_i = obj.sidx;
            sidx_j = agent2.sidx;
            
            %obj.sidx = (obj.sid-1)*obj.n_states+1;
            %obj.idxs = obj.sidx:(obj.sidx+obj.n_states-1);   
            
            x1 = obj.x(sidx_i);      y1 = obj.x(sidx_i+1); 
            x2 = obj.x(sidx_j);      y2 = obj.x(sidx_j+1); 
            %x2 = agent2.x(sidx_j);      y2 = agent2.x(sidx_j+1);
            r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );

            if r > 0
                H_line = zeros(1, obj.max_neighbors * obj.n_states);
                 
                H_line(sidx_i)        = (x1-x2)/r;
                H_line(sidx_i+1)      = (y1-y2)/r;
                H_line(sidx_j)        = (x2-x1)/r;
                H_line(sidx_j+1)      = (y2-y1)/r;

                obj.H_ext = [obj.H_ext; H_line]; 
                obj.z_ext = [obj.z_ext; d]; 

                obj.h_ext{length(obj.h_ext) + 1} = @(x) sqrt( ( x(sidx_i) - x(sidx_j) )^2 + ( x(sidx_i+1) - x(sidx_j+1) )^2 );
            end

          elseif strcmp(type, 'v2i')  
              
            xy_gt_i   = obj.gt(end, :);
            sidx_i = obj.sidx;  
            infra_node_xy = agent2(:, 2:3);
            obj.links = [obj.links; obj.id, agent2(1)]; 
            
            H_line = zeros(1, obj.max_neighbors * obj.n_states);        
            
            d = sqrt((xy_gt_i(1) - infra_node_xy(1))^2 + (xy_gt_i(2) - infra_node_xy(2))^2);
            d = d + normrnd(0, obj.system_setting.sigma_UWB);                  % Simulate ranging      
            
            x1 = obj.x(sidx_i);      y1 = obj.x(sidx_i+1); 
            x2 = infra_node_xy(1);   y2 = infra_node_xy(2); 
            r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );

            H_line(sidx_i)        = (x1-x2)/r;
            H_line(sidx_i+1)      = (y1-y2)/r;

            obj.H_ext = [obj.H_ext; H_line]; 
            obj.z_ext = [obj.z_ext; d]; 

            obj.h_ext{length(obj.h_ext) + 1} = @(x) sqrt( ( x(sidx_i) - x2 )^2 + ( x(sidx_i+1) - y2 )^2 );
            
          elseif strcmp(type, 'share-states')             

              %if isempty(find(obj.tracked == agent2.id, 1))
                 idxs_j = agent2.idxs;
                 obj.x(idxs_j) = agent2.x(idxs_j);
                 obj.P(idxs_j, idxs_j) = agent2.P(idxs_j, idxs_j);
                 obj.q(idxs_j) = agent2.q(idxs_j);
                 obj.Omega(idxs_j, idxs_j) = agent2.Omega(idxs_j, idxs_j);
                 obj.tracked = [obj.tracked; agent2.id];
                 %fprintf('Tracking: %i -> %i\n', obj.id, agent2.id);
              %end
              
          end
          
        end
        
        
        function reset(obj)
            obj.H_ext = []; 
            obj.z_ext = []; 
            obj.h_ext = {};
            obj.links = [];
        end
        
    end
end