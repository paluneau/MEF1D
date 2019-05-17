function [] = plot_mef1d(U,sol,lim_dom,err)
% PLOT_MEF1D
% Trace la sortie de la fonction SOLVE_MEF1D. Possibilit� de tracer la solution
% analytique pour comparer.
%
% !!! Param�tres !!!
%
% - U: la sortie de la fonctionSOLVE_MEF1D
%
% - sol (facultatif): la solution analytique du probl�me (function handle).
%
% - lim_dom (facultatif): domaine sur lequel on trace la solution
% analytique. Doit �tre obligatoirement sp�cifi� si sol est pass� en
% param�tre.
%
% - err (facultatif): Erreur en norme L^2 � afficher sur le graphique.

[~,m] = size(U);
figure;
hold on;
for i=1:m
    c1=plot(U{2,i},U{1,i}(U{2,i}),"r","LineWidth",1);
end

if nargin>1
    d = linspace(lim_dom(1),lim_dom(2),(lim_dom(2)-lim_dom(1))*30);
    c2=plot(d,sol(d),"b--","LineWidth",2);
    legend([c1 c2],"Num�rique (u_k)", "Analytique (u)");
end
if nargin==4
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',sprintf("||u-u_k||_{L^2} = %e",err),'Position',[0.15 0.8 0.1 0.1],'FitBoxToText','on','BackgroundColor','white');
end 
title("Solutions de p(x){\cdot}u(x)-d_x[q(x){\cdot}d_x[u(x)]]=R(x)");
xlabel("x");
ylabel("u(x)");
hold off;
end

