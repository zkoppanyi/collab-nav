function O = calc_observability(H, F)

        % Observaibility
        % See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.30.2399&rep=rep1&type=pdf
        O = H; 
        for l = 1 : size(H,1)
            O = [O; H*F^l]; 
        end
            