function [mu, C] = alpha_mtgp(logtheta,  x, y, xtest,x_train, f_train ,D,n_source)
% Predictions in Multi-task Gaussian process model
%
% [alpha, Kf, L, Kxstar, Kss] = alpha_mtgp(logtheta, covfunc_x, x, y, m, irank, ...
%					     nx, ind_kf, ind_kx, xstar )
%
% INPUT:
% - logtheta    : Vector of all parameters: [theta_lf; theta_x; sigma_l]
%                - theta_lf: the parameter vector of the
%                   cholesky decomposition of k_f
%                - theta_x: the parameters of K^x
%                - sigma_l: The log of the noise std deviations for each task
% - covfunc_x   : Name of covariance function on input space x
% - x           : Unique input points on all tasks 
% - y           : Vector of target values
% - m           : The number of tasks
% - irank       : The rank of K^f 
% - nx          : number of times each element of y has been observed 
%                usually nx(i)=1 unless the corresponding y is an average
% - ind_kx      : Vector containing the indexes of the data-points in x
%                which each observation y corresponds to
% - ind_kf      : Vector containing the indexes of the task to which
%                each observation y corresponds to
% - xstar       : input test points
% 
% OUTPUT:
%
% - alpha       : The solution to the ( (Kf x Kx) + (Sigma x I) )^{-1} y
% - Kf          : The Covariance matrix for the tasks
% - L           : The cholesky factorization of  C = (Kf x Kx) + Sigma
% - Kxstar      : Test-Train input covariances
% - Kss         : Test-test variances
% 
% Author: Edwin V. Bonilla
    
 % *** General settings here ****
%config = get_mtgp_config();
%MIN_NOISE = config.MIN_NOISE;
% ******************************

N = length(n_source);
for i = 1:N
    eval(['theta_lamda',num2str(i),'=[];'])
    eval(['theta_x',num2str(i),'=[];'])
    eval(['noise_source',num2str(i),'=[];'])
    eval(['noise_target',num2str(i),'=[];'])
    eval(['mu_t',num2str(i),'=[];'])
    eval(['C_t',num2str(i),'=[];'])
end

for i = 1:N
    eval(['theta_lamda',num2str(i),'=logtheta(i);'])
    for j = 1:D
        eval(['theta_x',num2str(i),'(j)=logtheta(N+j+(i-1)*D);'])
    end
    eval(['noise_source',num2str(i),'=logtheta(N+N*D+i);'])
    eval(['noise_target',num2str(i),'=logtheta(2*N+N*D+i);'])
end 
index =[1,n_source(1)];
for i =2:N
    index(i*2-1)=1+index(i*2-2);
    index(2*i)=index(i*2-2)+n_source(i);    
end

for i = 1:N
    Var_source = n_source(i);
    Var_target = size(xtest,1);
    source_x = x(index(2*i-1):index(2*i),:);
    target_x = xtest;
    source_y = y(index(2*i-1):index(2*i));
    %target_y = f_target_train;
    Kx11 = adptivecovSEard(eval(['theta_x',num2str(i)]), source_x);
    K11 = Kx11+(eval(['noise_source',num2str(i)])^2+1000*eps)*eye(Var_source);
    [Kx22, Kx12] = adptivecovSEard( eval(['theta_x',num2str(i)]), source_x, target_x);
    K22 = Kx22+(eval(['noise_target',num2str(i)])^2+eps)*eye(Var_target);
    K12 = Kx12*eval(['theta_lamda',num2str(i)]);
    K21 = K12';
    L =chol(K11)';
    alpha = solve_chol(L',source_y);
    eval(['mu_t',num2str(i),'=K21*alpha;'])
    v   = L\K21';
    eval(['C_t',num2str(i),'=K22-v''*v+eps*eye(Var_target);'])
   
end
wholeweight = 0;
for i= 1 :N
    eval(['wholeweight = wholeweight+abs(theta_lamda',num2str(i),');'])
end

mu=0;
for i =1:N
    eval(['mu=mu+abs(theta_lamda',num2str(i),')*mu_t',num2str(i),'/wholeweight;'])
end

C= 0;
for i =1:N
    eval(['C=C+(theta_lamda',num2str(i),')^2*sum(diag(C_t',num2str(i),'));'])
end


% ltheta_x = D+1;
% theta_lamda1 = logtheta(1);
% theta_x = zeros(D,1);
% for i = 1:(D)
% 
%       theta_x(i)=logtheta(i+1);
% end
% noise_source =logtheta(D+2);
% noise_target = logtheta(D+3);
% Var_source = size(x,1);
% Var_target = size(xtest,1);
% source_x = x;
% target_x = xtest;
% source_y = y;
% Kx11 = adptivecovSEard( theta_x, source_x);
% K11 = Kx11+(noise_source^2+100*eps)*eye(Var_source);
% [Kx22, Kx12] = adptivecovSEard( theta_x, source_x, target_x);
% K22 = Kx22+(noise_target^2+eps)*eye(Var_target);
% K12 = Kx12*theta_lamda1;
% K21 = K12';
% L =chol(K11)';
% alpha = solve_chol(L',source_y);
% mu_t = K21*alpha;
% v   = L\K21';
% C_t = K22-v'*v;
% 
% 
% 
% 
% 
% 


