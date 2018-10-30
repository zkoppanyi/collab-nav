function res = afun_eval(f, x)

    n = length(f);
    res = zeros(n,1);
    
    for i = 1 : n
        %fprintf('%i %i %i %i\n', i, length(res), length(f), length(x));
        res(i) = f{i}(x);
    end