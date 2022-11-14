% toy_example.m
% A toy example demonstrating how to use the mtgp package  for M=3 tasks 
%
% function toy_example()
% 
% 1. Generates sample from true MTGP model 
% 2. Assigns cell data for learning and prediction
% 3. Learns hyperparameters with minimize function
% 4. Makes predictions at all points on all tasks
% 5. Plots Data and predictions
%
% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
clear all; clc;
%rand('state',18);
%randn('state',20);
 
tic
PLOT_DATA       = 1;
%covfunc_x       = {'covSEard'};
M               = 5;    % Number of source domain
%D               = 10;    % Dimensionality of input space
% irank           = M;    % rank for Kf (1, ... M). irank=M -> Full rank


%% 1. Generating samples from true Model
%[x, Y, xtrain, ytrain, ind_kf_train, ind_kx_train , nx] = generate_data(covfunc_x, D, M);
%count = 1;
% target_k = [-1;-1.1917;-1.4281;-1.7320;-2.1445;-2.7475;-3.732;-5.6713;-11.4300;11.4300;5.6713;3.732;2.7475;2.1145;1.7320;1.4281;1.1917;1;0.8391;0.7002;0.5773;0.4663;0.3640;0.2680;0.1763;0.0875;-0.0875;-0.1763;-0.2680;-0.3640;-0.4663;-0.5773;-0.7002;-0.8391];
% angle = zeros(11,1);
% lamda1= zeros(11,1);
% lamda2  = zeros(11,1);
% xxx = zeros(11,1);
% map=zeros(11,1);
% getmin=zeros(11,1);
% getout=zeros(11,1);
% %nn =2;
% nn =length(target_k);
% la1 = zeros(nn,1);
% la2 = zeros(nn,1);
% ga1 = zeros(nn,1);
% ga2 = zeros(nn,1);
choose_mae = zeros(10,1);
choose_mse = zeros(10,1);
choose_mif = zeros(10,1);
for loop = 1: 10
x_multi_source=[];
f_multi_source=[];
n_multi_source=[];
for Nu = 1:M
[x_source,  f_source, x_target_train, x_target_test, f_target_train, f_target_test, D, n_source] = adaptivegenerate_data(Nu);
disp('generatedata_ok')
%target_k(jj)
x_multi_source = [x_multi_source;x_source];
f_multi_source = [f_multi_source;f_source];
n_multi_source = [n_multi_source;n_source];
end
%% 2. Assigns cell data for learning and prediction
%data  = {covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train, ind_kx_train};
data  = { x_multi_source, f_multi_source, x_target_train,  f_target_train, D,n_multi_source};

%% 3. Hyper-parameter learning
%[logtheta_all, deriv_range] = init_mtgp_default(M,D);
%[logtheta_all] = init_mtgp_default(M,D);
%disp('init_ok')
[logtheta,min_f,out_f]          = learn_mtgp( data);
disp('learn_ok')


%% 4. Making predictions at all points on all tasks
[ Ypred, Vpred ] = predict_mtgp_all_tasks(logtheta, data, x_target_test );
disp('predict_ok')

%% 5. Plotting the predictions
% if (PLOT_DATA)
%     plot_predictions(x_source, f_source,x_target_test, f_target_test , D, Ypred, Vpred);
%   %  figure, plot(nl); title('Negative Marginal log-likelihood');
% end
%disp(k)
%% 6. Normalized mean square error
[N,h] = size(Ypred);
sum =0;
sum1= 0;
for i =1 : N
    sum = sum+abs(f_target_test(i)-Ypred(i));
    sum1 = sum1+(f_target_test(i)-Ypred(i))^2;
end
MAE=sum/N;
MSE = sum1/N;
gap = f_target_test-Ypred;
%lamda =2*(1/(1+exp(-logtheta_all(3))))-1;
%lama = 2*(1/(1+logtheta(2)))^logtheta(1)-1;
%end
toc

choose_mae(loop) = MAE;
choose_mse(loop) = MSE;
choose_mif(loop) = min_f;
end
% tan = abs((k-1)/(1+k));
% at = atan(tan);
% angle(count)=at*180/pi;
% lamda1(count) = logtheta(1);
% lamda2 (count) = logtheta(2);
% xxx(count) = count;
% map(count)=sum;
% getmin(count)=min_f;
% getout(count)=out_f;
% count = count+1;
% %disp(k)
% disp(count)
% end
% [mg,inf]=min(getmin);
% la1(jj)=lamda1(inf);
% la2(jj)=lamda2(inf);
% ga1(jj)=mg;
% ga2(jj) = map(inf);
% % figure
% % plot3(lamda1,lamda2,angle)
% figure
% plot(xxx,map)
% end
% sc = [-90;-85;-80;-75;-70;-65;-60;-55;-50;-40;-35;-30;-25;-20;-15;-10;-5;0;5;10;15;20;25;30;35;40;50;55;60;65;70;75;80;85];
% %sc = [-90;-75];
% figure
% subplot(2,2,1);
% plot(sc,la1)
% subplot(2,2,2);
% plot(sc,la2)
% subplot(2,2,3);
% plot(sc,ga1)
% subplot(2,2,4);
% plot(sc,ga2)
% % 
% % 
% figure
% subplot(2,2,1);
% plot(target_k,la1)
% subplot(2,2,2);
% plot(target_k,la2)
% subplot(2,2,3);
% plot(target_k,ga1)
% subplot(2,2,4);
% plot(target_k,ga2)
