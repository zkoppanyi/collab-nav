function [x, P] = ukf_update( x, z, h, P, R, X1, X2)

        settings = ukf_settings();
        L           = numel(x);                                 %numer of states     
        m           = numel(z);                                 %numer of measurements
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
        
        [z1,Z1,P2,Z2] = ukf_unscented_transform(h,X1,Wm,Wc,m,R);       %unscented transformation of measurments
        P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
        K=P12*inv(P2);
        x=x+K*(z-z1);                              %state update
        P=P-K*P12';                                %covariance update
