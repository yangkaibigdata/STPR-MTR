function [c,cep]=fnonlon(t)
lamda1 = t(1);
lamda2 = t(2);
% theat1 =t(3);
% theat2 = t(4);
% noise_source = t(5);
% noise_target = t(6);
c = lamda1*lamda1-lamda2;
cep =[];
% if nargout>2
%     DC = [2*lamda1;-1;0;0;0;0];
%     DCep=[];
end
