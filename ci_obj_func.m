function res = ci_obj_func(betas, Is)

    n = length(betas);
    res = zeros(size(Is{1}, 1), size(Is{1}, 1));
    for i = 1 : n
        res = res + betas(i) * Is{i};
    end