function show_sparisty_pattern(J, fig_num, eps)

    if nargin < 2
        fig_num = 1;
    end

    if nargin < 3
        eps = -Inf;
    end

    sparI = J;
    sparI(sparI < eps) = 0;
    sparI(sparI == 0) = 0;
    sparI(sparI ~= 0) = 255;

    figure(fig_num); clf; hold on;
    image(flip(sparI));
    grid on;
    axis equal;