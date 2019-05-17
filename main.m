close all;clear all;
%% Probl�me type
% Trouver u(x) t.q. p(x)*u(x)-d_x[q(x)*d_x[u(x)]]=R(x) avec conditions aux
% limites de Dirichlet

%% Probl�me de laplacien
% Trouver u(x) t.q. -d^2_x[u(x)]=5 avec conditions aux
% limites de Dirichlet u(0)=0 et u(5)=2
p1=@(x)0;
q1=@(x)1;
R1=@(x)5; % En mettant une constante non-nulle, on s'assure d'avoir un polynome de degr� 2 comme solution
Omega1 = [0 5];
dirichlet1 = [0 2];

sol_analytique_laplacien = @(x)-5/2.*x.^2+129/10.*x;

[U1, err_L2] = solve_mef1d(p1,q1,R1,Omega1,1,dirichlet1,2,sol_analytique_laplacien);
plot_mef1d(U1,sol_analytique_laplacien,Omega1);

%% Probl�me de tension dans un c�ble
% Trouver u(x) t.q. -d_x[400*d_x[u(x)]]=r(x) avec conditions aux
% limites de Dirichlet u(0)=u(5)=0
p2=@(x)0;
q2=@(x)400; % Tension dans le cable
R2=@(x)r(x); % Poids lin�aire du cable + Poids ajout�, d�fini plus bas.
dirichlet2=[0 0];

[U2,~] = solve_mef1d(p2,q2,R2,Omega1,10,dirichlet2,2);
plot_mef1d(U2);

% Repr�sente le poids lin�aire du c�ble et une masse non constante r�partie
% sur la partie gauche du c�ble.
function y = r(x)
[~,dim] = size(x);
y = zeros(1,dim);
for i=1:dim
    if 0<=x(i) && x(i)<=2
        y(i)=-6-40*x(i);
    elseif 2<=x(i) && x(i)<=5
        y(i)=-6;
    end
end
end