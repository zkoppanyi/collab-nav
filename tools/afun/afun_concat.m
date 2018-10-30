function f = afun_concat(f1, f2)

    f= {};
    for i = 1 : length(f1)
        f{i} = f1{i};
    end
    
    for i = 1 : length(f2)
        f{length(f1)+i} = f2{i};
    end
    