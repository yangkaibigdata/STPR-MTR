function [x_source, f_source, x_target_train, x_target_test, f_target_train, f_target_test, D,n_source] = adaptivegenerate_data(Nu)
% [x, Y, xtrain, ytrain, ind_kx_train, ind_kf_train, nx] = ...
%                                    generate_data(covfunc_x, D, M)
% Generates data for testing MTGP model
%
% INPUT:
% - covfunc_x : Input covariance function
% - D         : Input dimensionality
% - M         : Number of tassks
%
% OUTPUT:    
% - x             : all inputs
% - Y             : all target values
% - xtrain        : Training inputs
% - ytrain        : Training target values
% - ind_kf_train  : Vector containing the indexes of the task to which
%                each observation y corresponds to
% - ind_kx_train  : Vector containing the indexes of the data-points in x
%                which each observation y corresponds to
% - nx            : Number of observations per each task input
%
% Edwin V. Bonilla


% rho      = 0.8;    % Correlations between problems
% sigma_2n = 0.0001; % Noise variance
% 
% 
% %% Input Point
% if (D==1)
%     range = linspace(-1,1, 200)'; % Alternative
%     x = range;
% elseif (D==2)    
%    range = linspace(-1,1, 20)'; % Alternative
%   [X1,X2] = meshgrid(range, range);%[A B]=meshgrid(a b),AµÈÓÚa
%   x = [X1(:), X2(:)];
% end
% [N D] = size(x);    % Total number of samples per task
% n = N*M;
% 
% 
% %% Task covariance
% Kf = rho*ones(M);
% idx_diag = sub2ind(size(Kf), 1:M, 1:M);%´ÓÁÐ¿ªÊ¼Ë÷Òý£¬½«kf×ª»¯³É¶Ô½Ç¾ØÕó
% Kf(idx_diag) = 1;
% 
% %% Input covariance
% ltheta_x =  eval(feval(covfunc_x{:}));
% theta_x = log(ones(ltheta_x,1));      
% 
% %% Noise variances (re-parametrized)
% theta_sigma = log(sqrt(sigma_2n*ones(M,1))); 
% 
% %% Full covariance
% Sigma2 = diag(exp(2*theta_sigma));
% Kx = feval(covfunc_x{:}, theta_x, x);
% C = kron(Kf,Kx) + kron(Sigma2,eye(N))
% %C = kron(Kf,Kx)
% L = chol(C)';    
% 
% %% Observations
% y      = L*randn(n,1); % Noisy observations
% v      = repmat((1:M),N,1); 
% ind_kf = v(:);  % indices to task
% v      = repmat((1:N)',1,M);
% ind_kx = v(:);           % indices to input space data-points
% Y      = reshape(y,N,M); % Matrix of observations across all tasks
% 
% %% Selecting data-points for training
% ntrain       = floor(0.1*n);
% v            = randperm(n);
% idx_train    = v(1:ntrain);
% nx           = ones(ntrain,1); % observations on each task-input point
% ytrain       = y(idx_train);
% xtrain       = x; 
% ind_kx_train = ind_kx(idx_train);
% ind_kf_train = ind_kf(idx_train);

%% target domain
% MU=0;
% SIGMA=0.2;
% m=200;
% n=1;
% R = normrnd(MU,SIGMA,m,n);
% %w0=1+100*rand(100,1);
% rand('state',20);
% D=1;
% num=10;
% %w0=randperm(num,D)';
% %w0 =[2;3];
% % %w0=[17;19;70;50;95;11;93;90;36;55;10;40;89;85;97;38;35;52;84;26;6;94;79;72;53;43;33;78;21;12;80;60;9;99;92;2;56;27;76;30;16;88;69;65;37;22;14;54;100;58;13;81;68;47;44;62;75;59;64;23;42;4;45;87;29;77;73;18;15;41;49;61;83;48;66;74;28;63;3;46;34;96;67;91;1;32;98;71;24;20;5;31;8;51;86;25;57;82;39;7];
% a=1;
% b=10;
% row=50;
% % for i = 1:100
% %     x_target(i)=i;
% % end
% % x_target = x_target';
% range = linspace(-40,40, 200)';
% %range2 = linspace(30,60, 200)';
% x_target = range;
% % % x_target(:,1) = range;
% % % x_target(:,2) = range2;
% % %x_target=a+b*rand(row,D);
% % %f_target=x_target*w0
% % %Z = R(151:300);
% % %size(x_target);
% f_target=x_target*2+R;
% %f_target=R;
% % 
% % range = linspace(-1,1, 100)'; % Alternative
% % x_target = range;
% % [row, col]=size(x_target);
% % D= col+1;
% % theta_x = [log(2);log(3)]; 
% % Kx = adptivecovSEard( theta_x, x_target);
% % sigma_noise = 0.3;
% % C = Kx+sigma_noise*eye(row);
% % L = chol(C)';
% % s = randn(row,1);
% % f_target = L*s;
% 
%  %x_target =linspace(-2,2, 40)';
% % f_target = (x_target-1).*(x_target-1);
% pick=10;
% % D=1;
%  [row,col] = size(x_target);
% %index_f=[30 50 70 90 110 130 150 170 200 240];
%  index_f=randperm(row,pick);
% x_target_train=x_target(index_f,:);
% x_target_test=x_target;
% x_target_test(index_f,:)=[];
% 
% f_target_train=f_target(index_f);
% f_target_test=f_target;
% f_target_test(index_f)=[];


% target=load('wine-red.txt');
% x_target=target(1:20,1:11);
% f_target=target(1:20,12);
% [a, D]= size(x_target);
% pick=floor(a*0.05);
% index_f=randperm(a,pick);
% x_target_train=x_target(index_f,:);
% x_target_test=x_target;
% x_target_test(index_f,:)=[];
% f_target_train=f_target(index_f);
% f_target_test=f_target;
% f_target_test(index_f)=[];

%% source domain
% delta=0;
% a1=15;
% b1=30;
% %x_source=a1+b1*rand(row,D);
% 
% f_source=x_target*w0+R;
%w=a1+b1*rand(D,1);
%g=x*(w0+delta*w)+R;
%x_source = x_target;
%f_source = f_target;
%f_source =x_source*(w0+2.199);
%f_source =x_source*8+R;
% x_source = linspace(-3,3, 100)';
%f_source = f_target;
% w1 =[4;5];
 %f_source =x_source*3+R;
 %f_source = sin(x_source).*(x_source).*(x_source)+R;
%f_source = (x_source+1).*(x_source+1).*(x_source+1)+R;
 % f_source = -5*(x_source-1).*(x_source-1)+R+1;
% x_source = x_target;
% theta_source = [log(2.5);log(3)]; 
% Kxx = adptivecovSEard( theta_source, x_source);
% sigma_noisesource = 0.4;
% C1 = Kxx+sigma_noisesource*eye(row);
% L1 = chol(C1)';
% f_source = L1*s;

% source=load('wine-white.txt');
% x_source=source(1:100,1:11);
% f_source=source(1:100,12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D=1;
% x_source = linspace(0,10,100)';
% f_source =10*sin(x_source);
% x_target_test = linspace(0,10,40)';
% f_target_test =2*x_target_test-8;
% pick =[1,5,10,15,20,25,30,35,40];
% x_target_train = x_target_test(pick);
% f_target_train = f_target_test(pick);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D=2;
% t=-2:0.1:2;
% [x,y] =meshgrid(t);
% 
% a = reshape(x,1,1681);
% x_source(1,:) = a;
% b = reshape(y,1,1681);
% x_source(2,:) = b;
% f_source =10*( x_source(1,:)+x_source(2,:).^5+x_source(1,:).^3).*exp(-x_source(1,:).^2-x_source(2,:).^2);
% 
% x_target_test(1,:) = a;
% x_target_test(2,:) = b;
% f_target_test = -8*(x_target_test(1,:)+x_target_test(1,:).^3).*exp(-x_target_test(1,:).^2 - x_target_test(2,:).^2);
% num =1;
% for j = 1:13
%     for i = 1:13
% pick(num) = j*(3*41)+i*3;
% num = num+1;
%     end
% end
% 
% x_target_train =x_target_test(:,pick);
% f_target_train =f_target_test(pick);
% 
% x_source =x_source';
% f_source = f_source';
% x_target_train=x_target_train';
% f_target_train = f_target_train';
% x_target_test = x_target_test';
% f_target_test = f_target_test';
% n_source = 1681;
% n_target_train = num-1;
% n_target_test = 1681;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Nu ==1
% x_source = linspace(-10,10,100);
% f_source =-x_source.^2+20;
% x_target_test = linspace(-10,10,40);
% f_target_test = 3*x_target_test.^2;
% pick =[1,5,10,15,20,25,30,35,40];
% x_target_train = x_target_test(pick);
% f_target_train = f_target_test(pick);
%  [D,n_source] = size(x_source);
% % [~,n_target_train] = size(x_target_train);
% % [~,n_target_test] = size(x_target_test);
% x_source =x_source';
% f_source = f_source';
% x_target_train=x_target_train';
% f_target_train = f_target_train';
% x_target_test = x_target_test';
% f_target_test = f_target_test';
% end
% if Nu ==2
% x_source = linspace(-10,10,100);
% f_source =x_source.^2-20;
% x_target_test = linspace(-10,10,40);
% f_target_test = 3*x_target_test.^2;
% pick =[1,5,10,15,20,25,30,35,40];
% x_target_train = x_target_test(pick);
% f_target_train = f_target_test(pick);
% [D,n_source] = size(x_source);
% % [~,n_target_train] = size(x_target_train);
% % [~,n_target_test] = size(x_target_test);
% x_source =x_source';
% f_source = f_source';
% x_target_train=x_target_train';
% f_target_train = f_target_train';
% x_target_test = x_target_test';
% f_target_test = f_target_test';
% 
% end
% if Nu ==3
% x_source = linspace(-10,10,100);
% f_source =2*x_source.^2+x_source-10;
% x_target_test = linspace(-10,10,40);
% f_target_test = 3*x_target_test.^2;
% pick =[1,5,10,15,20,25,30,35,40];
% x_target_train = x_target_test(pick);
% f_target_train = f_target_test(pick);
% [D,n_source] = size(x_source);
% % [~,n_target_train] = size(x_target_train);
% % [~,n_target_test] = size(x_target_test);
% x_source =x_source';
% f_source = f_source';
% x_target_train=x_target_train';
% f_target_train = f_target_train';
% x_target_test = x_target_test';
% f_target_test = f_target_test';
% 
% end
% if Nu ==4
% x_source = linspace(-10,10,100);
% f_source =0.5*x_source.^2-x_source-10;
% x_target_test = linspace(-10,10,40);
% f_target_test = 3*x_target_test.^2;
% pick =[1,5,10,15,20,25,30,35,40];
% x_target_train = x_target_test(pick);
% f_target_train = f_target_test(pick);
% [D,n_source] = size(x_source);
% % [~,n_target_train] = size(x_target_train);
% % [~,n_target_test] = size(x_target_test);
% x_source =x_source';
% f_source = f_source';
% x_target_train=x_target_train';
% f_target_train = f_target_train';
% x_target_test = x_target_test';
% f_target_test = f_target_test';
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traindata = csvread('deal_train_scarcos.csv');
testdata = csvread('deal_train_scarcos.csv');
source_1 = [1,8,15];
source_2 = [2,9,16];
source_3 = [3,10,17];
source_4 = [4,11,18];
source_5 = [5,12,19];
source_6 = [6,13,20];
source_7 = [7,14,21];
target_1 = 22;
target_2 = 23;
target_3 = 24;
target_4 = 25;
target_5 = 26;
target_6 = 27;
target_7 = 28;
a =[1:100:2500];
% t = randperm(2500,50);
% s = randperm(12500,500);
t = 1:500;

s = 25:50;
if Nu ==1
    x_target_train = traindata(s,source_1);
    %y_target_train = traindata(25:50,522);
    f_target_train = traindata(s,target_1);
    [a, D]= size(x_target_train);
    x_source = traindata(t,source_2);
    f_source = traindata(t,target_2);
    %n_target_train = length(f_target_train);
    n_source =length(f_source);
    x_target_test = testdata(2000:2999,source_1);
    f_target_test = testdata(2000:2999,target_1);
    %n_target_test = length(f_target_test);
end

if Nu ==2
    x_target_train = traindata(s,source_1);
    %y_target_train = traindata(25:50,522);
    f_target_train = traindata(s,target_1);
    [a, D]= size(x_target_train);
    x_source = traindata(t,source_4);
    f_source = traindata(t,target_4);
    %n_target_train = length(f_target_train);
    n_source =length(f_source);
    x_target_test = testdata(2000:2999,source_1);
    f_target_test = testdata(2000:2999,target_1);
    %n_target_test = length(f_target_test);
end

if Nu ==3
    x_target_train = traindata(s,source_1);
    %y_target_train = traindata(25:50,522);
    f_target_train = traindata(s,target_1);
    [a, D]= size(x_target_train);
    x_source = traindata(t,source_5);
    f_source = traindata(t,target_5);
    %n_target_train = length(f_target_train);
    n_source =length(f_source);
    x_target_test = testdata(2000:2999,source_1);
    f_target_test = testdata(2000:2999,target_1);
    %n_target_test = length(f_target_test);
end

if Nu ==4
    x_target_train = traindata(s,source_1);
    %y_target_train = traindata(25:50,522);
    f_target_train = traindata(s,target_1);
    [a, D]= size(x_target_train);
    x_source = traindata(t,source_6);
    f_source = traindata(t,target_6);
    %n_target_train = length(f_target_train);
    n_source =length(f_source);
    x_target_test = testdata(2000:2999,source_1);
    f_target_test = testdata(2000:2999,target_1);
    %n_target_test = length(f_target_test);
end

if Nu ==5
    x_target_train = traindata(s,source_1);
    %y_target_train = traindata(25:50,522);
    f_target_train = traindata(s,target_1);
    [a, D]= size(x_target_train);
    x_source = traindata(t,source_7);
    f_source = traindata(t,target_7);
    %n_target_train = length(f_target_train);
    n_source =length(f_source);
    x_target_test = testdata(2000:2999,source_1);
    f_target_test = testdata(2000:2999,target_1);
    %n_target_test = length(f_target_test);
end


return;



