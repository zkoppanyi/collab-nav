function f = afun_create(n, f_proto)

    f= {};
    for i = 1 : n
        f{i} = f_proto;
    end