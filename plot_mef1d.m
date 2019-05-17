function [] = plot_mef1d(U,sol,lim_dom)
% PLOT_MEF1D
% Trace la sortie de la fonction SOLVE_MEF1D. Possibilité de tracer la solution
% analytique pour comparer.
%
% !!! Paramètres !!!
%
% - U: la sortie de la fonctionSOLVE_MEF1D
%
% - sol (facultatif): la solution analytique du problème (function handle).
%
% - lim_dom (facultatif): domaine sur lequel on trace la solution
% analytique. Doit être obligatoirement spécifié si sol est passé en
% paramètre.

[~,m] = size(U);
figure;
hold on;
for i=1:m
    c1=plot(U{2,i},U{1,i}(U{2,i}),"r","LineWidth",1);
end

if nargin>1
    d = linspace(lim_dom(1),lim_dom(2),(lim_dom(2)-lim_dom(1))*30);
    c2=plot(d,sol(d),"b--","LineWidth",2);
    legend([c1 c2],"Numérique", "Analytique");
end
title("Solutions de p(x)*u(x)-d_x[q(x)*d_x[u(x)]]=R(x)");
xlabel("x");
ylabel("u(x)");
hold off;
end

