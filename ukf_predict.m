function [x, P, X1, X2] = ukf_predict(f, x, P, Q)

        settings = ukf_settings();
        
        L           = numel(x);                                 %numer of states       
        alpha       = settings.alpha;                           %default, tunable
        ki          = settings.ki;                              %default, tunable
        beta        = settings.beta;                            %default, tunable
        lambda      = alpha^2*(L+ki)-L;                         %scaling factor
        c           = L+lambda;                                 %scaling factor
        Wm          = [lambda/c 0.5/c+zeros(1,2*L)];            %weights for means
        Wc          = Wm;
        Wc(1)       = Wc(1)+(1-alpha^2+beta);                   %weights for covariance
        c           = sqrt(c);
        X           = ukf_sigma_points(x,P,c);                  %sigma points around x
        [x, X1, P, X2] = ukf_unscented_transform(f,X,Wm,Wc,L,Q);   
