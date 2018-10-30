clear variables;
clc; close all;

P11 = [0.7 0.4; 0.4 0.3];
P22 = [0.2 -0.1; -0.1 0.4];

figure(1); clf; hold on;
error_ellipse(P11, 'style', 'k');
error_ellipse(P22, 'style', 'k');
for i = 1 : 10
    %P12 = [0 0; 0 0]/3;
    %P12 = (rand(2)-0.5)+(rand(2)-0.5);
    P12 = (rand(2)-0.5)/6;
    Pf = inv( inv(P11) + (inv(P11)*P12 - eye(2)) * inv(P22 - P12'*inv(P11)*P12) * (P12'*inv(P11) - eye(2)) );

    P = [P11, P12; P12', P22];
    if min(eig(P)) < 0
        continue
    end    
    error_ellipse(Pf,  'style', 'b-');
    axis equal
end

w0 = [0.1 0.9];
Aeq     = [1 1];
beq     = 1;
lb      = [0 0];
ub      = [1 1];    
Pci = @(w) inv(w(1)*inv(P11) + w(2)*inv(P22));
w = fmincon(@(w) trace(Pci(w)), w0, [], [], Aeq, beq, lb, ub);         
error_ellipse(Pci(w), 'style', 'r');
error_ellipse(Pci([0.5 0.5]), 'style', 'g');
grid on;
set(gca, 'FontSize', 14)

