classdef SelfishAgent < handle
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
        R_ext    = [];
        Q        = [];
        
        links    = [];
        
        system_setting = {};
        
    end
    
    methods        
        
        function obj = SelfishAgent(id, x_init, P_init, t_init, max_neighbors, system_setting)
            
            persistent nsidx
            if isempty(nsidx)
               nsidx = 1;
            else
                nsidx = nsidx + 1;
            end
            
            obj.id              = id;
            obj.sid             = nsidx;
            obj.max_neighbors   = max_neighbors;
            obj.n_states        = length(x_init);
            obj.x               = x_init;
            obj.P               = P_init;
            obj.t0              = t_init;
            obj.system_setting  = system_setting;
            
        end
        
        function apply_update(obj, x, P)
            
                if (size(x, 1) > size(x, 2))
                    x2 = x';
                end
            
                obj.x = x;
                obj.P = P;
                obj.x_hist = [obj.x_hist; x2];
        end
            
        function build_predict(obj)

            obj.F = [1 0 sin(obj.x(4))*obj.system_setting.dt 0; 0 1 cos(obj.x(4))*obj.system_setting.dt 0; 0 0 1 0; 0 0 0 1];

            obj.f = cell(4, 1);
            obj.f{1} = @(x) x(1) +  x(3)*sin(x(4)) * obj.system_setting.dt;
            obj.f{2} = @(x) x(2) +  x(3)*cos(x(4)) * obj.system_setting.dt;
            obj.f{3} = @(x) x(3);
            obj.f{4} = @(x) x(4);

            obj.Q = diag([0.5, 0.5, 0.05, (2.5/180*pi)^2]);
        end
        
        function build_int_update(obj, type, upt)
            
            if (size(upt, 1) > size(upt, 2))
                upt = upt';
            end
            
            obj.gt = [obj.gt; upt];
        
            if or(or( strcmp(type, 'GPS'), strcmp(type, 'IMU')), strcmp(type, 'GPS+IMU'))       
                H = zeros(obj.n_states, obj.n_states);
                R = eye(obj.n_states, obj.n_states);
                z = zeros(obj.n_states,1);
                h = cell(obj.n_states, 1);
            else
                H = obj.H;
                z = obj.z;
            end

            if strcmp(type, 'GPS')

                H(1,1) = 1; 
                H(2,2) = 1; 
                
                R(1,1) = obj.system_setting.sigma_GPS^2; 
                R(2,2) = obj.system_setting.sigma_GPS^2;

                h{1}   =  @(x) x(1);
                h{2}   =  @(x) x(2);        
                h{3}   =  @(x) 0;
                h{4}   =  @(x) 0;

                z(1)     =  upt(1) + normrnd(0, obj.system_setting.sigma_GPS);        
                z(2)     =  upt(2) + normrnd(0, obj.system_setting.sigma_GPS);                 
                z(3)     =  0;        
                z(4)     =  0; 

            elseif strcmp(type, 'IMU')

                H(3,3) = 1; 
                H(4,4) = 1; 
                R(3,3) = obj.system_setting.sigma_v^2; 
                R(4,4) = obj.system_setting.sigma_IMU^2;

                h{1}     =  @(x) 0;
                h{2}     =  @(x) 0;
                h{3}     =  @(x) x(3);
                h{4}     =  @(x) x(4);

                z(1)     =  0;        
                z(2)     =  0; 
                z(3)     =  upt(3) + normrnd(0, obj.system_setting.sigma_v); 
                z(4)     =  upt(4) + normrnd(0, obj.system_setting.sigma_IMU);  

            elseif strcmp(type, 'GPS+IMU')

                H = eye(4);
                R(1,1) = obj.system_setting.sigma_GPS^2; 
                R(2,2) = obj.system_setting.sigma_GPS^2;
                R(3,3) = obj.system_setting.sigma_v^2; 
                R(4,4) = obj.system_setting.sigma_IMU^2;

                h{1}     =  @(x) x(1);
                h{2}     =  @(x) x(2);
                h{3}     =  @(x) x(3);
                h{4}     =  @(x) x(4);

                z(1)     =   upt(1) + normrnd(0, obj.system_setting.sigma_GPS);        
                z(2)     =   upt(2) + normrnd(0, obj.system_setting.sigma_GPS);  
                z(3)     =   upt(3) + normrnd(0, obj.system_setting.sigma_v); 
                z(4)     =   upt(4) + normrnd(0, obj.system_setting.sigma_IMU);           

            end

            obj.H_int = H;
            obj.h_int = h;
            obj.z_int = z;
            obj.R_int = R;
            
            % Save some history data             
            obj.z_int_hist = [obj.z_int_hist; z'];     
            
        end
        
        function comm(obj, agent2, type, varargin)            
          
          if strcmp(type, 'UWB') 
            
            xy_gt   = obj.gt(end, :);
            xy_gt_j = agent2.gt(end, :);
            
            % This means that the vehicle left AOI, so the coordinates are
            % not updated; skip this communication link
            if ~(and(and(and(xy_gt_j(1) > obj.system_setting.AOI(1,1), ...
                    xy_gt_j(1) < obj.system_setting.AOI(1,2)), ...
                    xy_gt_j(2) > obj.system_setting.AOI(2,1)), ...
                    xy_gt_j(2) < obj.system_setting.AOI(2,2)))
                
                return;
            end
            
            obj.links = [obj.links; obj.id, agent2.id]; 
            d = sqrt((xy_gt(1) - xy_gt_j(1))^2 + (xy_gt(2) - xy_gt_j(2))^2);
            d = d + normrnd(0, obj.system_setting.sigma_UWB);                  % Simulate ranging            
            
            H_line = zeros(1, obj.max_neighbors * obj.n_states);
        
            sidx = (obj.sid - 1)*obj.n_states+1;
            sidx_j = (agent2.sid - 1)*obj.n_states+1;

            x1 = obj.x(1);      y1 = obj.x(2); 
            x2 = agent2.x(1);   y2 = agent2.x(2); 
            r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );

            H_line(sidx)        = (x1-x2)/r;
            H_line(sidx+1)      = (y1-y2)/r;
            H_line(sidx_j)      = (x2-x1)/r;
            H_line(sidx_j+1)    = (y2-y1)/r;

            obj.H_ext = [obj.H_ext; H_line]; 
            obj.z_ext = [obj.z_ext; d]; 
            %H_x = [H_x; r];          %#ok

            obj.h_ext{length(obj.h_ext) + 1} = @(x) sqrt( ( x(sidx) - x(sidx_j) )^2 + ( x(sidx+1) - x(sidx_j+1) )^2 );
            
        elseif strcmp(type, 'v2i')  
              
            xy_gt_i   = obj.gt(end, :);
            sidx_i = (obj.sid - 1)*obj.n_states+1;  
            infra_node_xy = agent2(:, 2:3);
            obj.links = [obj.links; obj.id, agent2(1)]; 
            
            H_line = zeros(1, obj.max_neighbors * obj.n_states);        
            
            d = sqrt((xy_gt_i(1) - infra_node_xy(1))^2 + (xy_gt_i(2) - infra_node_xy(2))^2);
            d = d + normrnd(0, obj.system_setting.sigma_UWB);                  % Simulate ranging      
            
            x1 = obj.x(1);           y1 = obj.x(2); 
            x2 = infra_node_xy(1);   y2 = infra_node_xy(2); 
            r = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );

            H_line(sidx_i)        = (x1-x2) / r;
            H_line(sidx_i+1)      = (y1-y2) / r;

            obj.H_ext = [obj.H_ext; H_line]; 
            obj.z_ext = [obj.z_ext; d]; 

            obj.h_ext{length(obj.h_ext) + 1} = @(x) sqrt( ( x(sidx_i) - x2 )^2 + ( x(sidx_i+1) - y2 )^2 );

          elseif strcmp(type, 'share-states')             

              if isempty(find(obj.tracked == agent2.id, 1))
                 idxs_j = agent2.idxs;
                 obj.x(idxs_j) = agent2.x(idxs_j);
                 obj.P(idxs_j, idxs_j) = agent2.P(idxs_j, idxs_j);
                 obj.tracked = [obj.tracked; agent2.id];
                 %fprintf('Tracking: %i -> %i\n', obj.id, agent2.id);
              end
                
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