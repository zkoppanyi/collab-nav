function dr = obj_fn(x, x0l, meas, ci)

        r = @(x1, y1, x2, y2) (x1-x2).^2 + (y1-y2).^2;
        dr = nan(size(meas, 1), 1);
        for k = 1 : size(meas, 1)
            
            idx_i   = meas(k,1);
            idx_j   = meas(k,2);
            from_i  = (idx_i-1)*2 + 1;
            from_j  = (idx_j-1)*2 + 1;
            
            p1 = isempty(find(ci == idx_i, 1));
            p2 = isempty(find(ci == idx_j, 1));
            
            %[from_i from_i+1 from_j from_j+1 sqrt(r( x(from_i), x(from_i+1), x(from_j), x(from_j+1) )) meas(k,3)]

            if and(~p1, p2)
                dr(k) = sqrt(r( x0l(from_i), x0l(from_i+1), x(from_j), x(from_j+1) )) - meas(k,3);
            elseif and(p1, ~p2)
                dr(k) = sqrt(r( x(from_i), x(from_i+1), x0l(from_j), x0l(from_j+1) )) - meas(k,3);
            elseif and(~p1, ~p2)               
                dr(k) = 0;
            else
                dr(k) = sqrt(r( x(from_i), x(from_i+1), x(from_j), x(from_j+1) )) - meas(k,3);
            end
        end
        