
% This function can display conditional data in simulation grid.

% The cnd_data can be a 2-D matrix with arbitary size for example 250*250. 
% ---------------- Programer : Mojtaba Bavandsavadkoohi------------------

function [zero_element, non_zero_element] = Display_cnd_data(cnd_data)
[m, n] = size(cnd_data);
zero_element = find(cnd_data==0);
non_zero_element = find(cnd_data==1);
[I1, J1] = ind2sub([m n],zero_element);
[I2, J2] = ind2sub([n n],non_zero_element);
axis square 
% scatter(J1,I1,'o' , 'MarkerEdgeColor','c')
% scatter(J1,I1,'filled','ks','MarkerEdgeColor','r','LineWidth',.1)
scatter(J1,I1,'o' , 'MarkerEdgeColor','c','LineWidth',1);
hold on
axis square 
xlim([0 m]);ylim([0 n]);xticks(0:200:m);yticks(0:200:n)
set(gca,'FontSize',12);
% scatter(J2,I2,'o' , 'MarkerEdgeColor','c')
% scatter(J2,I2,'filled','ws','MarkerEdgeColor','b','LineWidth',.1)
scatter(J2,I2,'o' , 'MarkerEdgeColor','r ','LineWidth',1);
%r
