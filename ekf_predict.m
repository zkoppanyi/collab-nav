function [x_hat, P_hat] = ekf_predict(f, F, x, P, Q)
        
        if ~isempty(f)
            if iscell(f)
                x_hat = afun_eval(f, x);
            else
                x_hat = f(x);
            end            
        else
             x_hat = F*x;
        end
        
        P_hat = F*P*F' + Q; 