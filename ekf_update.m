function [x, P, K] = ekf_update( x, z, h, H, P, R)

        if ~isempty(h)
            if iscell(h)
                y = z - afun_eval(h, x);
            else
                y = z - h(x);
            end            
        else
             y = z - H*x;
        end
        
        S = R + H*P*H';
        K = P*H'*inv(S);
        x = x + K*y;
        P = (eye(size(K, 1), size(K, 1)) - K*H)*P*(eye(size(K, 1), size(K, 1)) - K*H)'+K*R*K';
        y = z - H*x;
        