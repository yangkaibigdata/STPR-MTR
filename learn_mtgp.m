function [logtheta,min_f,out_f] = learn_mtgp( data )
%LEARN_MTGP Learns hyperparameters of mtgp model using minimize
%
% INPUT:
% - logtheta_all : All initial hyper-parameter values 
% - deriv_range  : Indices of hyper-parameters to optimize for
% - data         : cell data in the order 
%                  [covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train, ind_kx_train]
%
% OUTPUT:
% - logtheta_all : all learned hyperparameters
% - nl           : Final negative marginal likelihood
%
% Edwin V. Bonilla

%% This can be changed: See minimize function for details
%niter     = 1000; % setting for minimize function: number of function evaluations


% ************* Learning Hyperparameters here *******************************
%logtheta0 = logtheta_all(deriv_range); % selects hyper-parameters to optimize
[ x_source, f_source, x_target_train, f_target_train, D,n_source] = deal(data{:});
% disp(D)
% x0 = init_mtgp_default(D);
x0 = init_mtgp_default(D,n_source);
% disp('init_ok')
%[0;0.1;0.01;0.01;0.01;(1e-7)*rand(2,1)];
%x0 = [0;0.01;0.01;(1e-7)*rand(2,1)];
%options= optimoptions('fmincon','display','iter');
%options= optimoptions('fmincon','Algorithm','active-set','GradObj','on','GradConstr','on','display','iter');
%[x,fval,exitflag,output] = fmincon(@(t) object_function(t,x_source,f_source,x_target_train,f_target_train,D),x0,[],[],[],[],[-1],[1],@(t) fnonlon(t));
[x,fval,exitflag,output] = fmincon(@(t) object_function(t,x_source,f_source,x_target_train,f_target_train,D,n_source),x0,[],[],[],[],[-1],[1]);
%[x,fval,exitflag,output] = fmincon(@(t) object_function(t,x_source,f_source,x_target_train,f_target_train),x0,options);
%[logtheta, nl] = minimize(logtheta0,'nmargl_mtgp',niter, logtheta_all, ...
%			 covfunc_x, x_source, f_source, x_target_train, f_target_train, D,deriv_range);

%% Update whole vector of parameters with learned ones
%logtheta_all(deriv_range) = logtheta;
logtheta = x;
min_f = fval;
out_f = exitflag;


return;



