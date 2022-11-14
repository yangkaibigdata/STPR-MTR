function [f]=object_function(t,x_source,f_source,x_target_train,f_target_train,D,n_source)
%theta_lamda1 = 0.96;
%disp('get in')
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
    eval(['theta_lamda',num2str(i),'=t(i);'])
    for j = 1:D
        eval(['theta_x',num2str(i),'(j)=t(N+j+(i-1)*D);'])
    end
    eval(['noise_source',num2str(i),'=t(N+N*D+i);'])
    eval(['noise_target',num2str(i),'=t(2*N+N*D+i);'])
end

% theta_lamda1 = t(1);
% theta_x = zeros(D,1);
% for i = 1:(D)
%     %theta_x(i)=log(abs(t(i+2)))
%       theta_x(i)=t(i+1);
% end
% noise_source =t(D+2);
% noise_target = t(D+3);
index =[1,n_source(1)];
for i =2:N
    index(i*2-1)=1+index(i*2-2);
    index(2*i)=index(i*2-2)+n_source(i);    
end

for i = 1:N
    Var_source = n_source(i);
    Var_target = size(x_target_train,1);
    source_x = x_source(index(2*i-1):index(2*i),:);
    target_x = x_target_train;
    source_y = f_source(index(2*i-1):index(2*i));
    target_y = f_target_train;
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
% Var_source = size(x_source,1);
% Var_target = size(x_target_train,1);
% source_x = x_source;
% target_x = x_target_train;
% source_y = f_source;
% target_y = f_target_train;
%if theta_lamda1^2<=theta_lamda2
% if theta_lamda1^2>theta_lamda2
%      theta_lamda2 = theta_lamda1^2+0.1;
%  end
%ep = 1;
%size(source_x)
% Kx11 = adptivecovSEard(theta_x, source_x);
% %L =chol(Kx11)'
% K11 = Kx11+(noise_source^2+1000*eps)*eye(Var_source);
% %K11 = Kx11+ep*eye(Var_source);
% [Kx22, Kx12] = adptivecovSEard( theta_x, source_x, target_x);
% K22 = Kx22+(noise_target^2+eps)*eye(Var_target);
% %K22 = Kx22+ep*eye(Var_target);
% %lamda = 2*(1/(1+lamda_mu))^lamda_b-1;
% K12 = Kx12*theta_lamda1;
% %K12 = Kx12*(2*(1/(1+theta_mu))^theta_b-1);
% K21 = K12';
% 
% 
% L =chol(K11)';
% alpha = solve_chol(L',source_y);
% %mu_target = sum(target_y)/Var_target;
% mu_t = K21*alpha;
% v   = L\K21';
% 
% C_t = K22-v'*v+eps*eye(Var_target);
% 
% Lt = chol(C_t)';
% y = target_y-mu_t;
% alpha_t = solve_chol(Lt',y);
% 
% n=Var_source+Var_target;
% 
% % negative log-likelihood
% f = 0.5*y'*alpha_t + sum(log(diag(Lt))) + 0.5*n*log(2*pi);
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
    eval(['C=C+((theta_lamda',num2str(i),')/wholeweight)^2*(C_t',num2str(i),');'])
end

f=sum(((target_y-mu)+diag(C)).^2);
%f=sum((target_y-mu).^2);
% else 
%     f = Inf;
% end
% if nargout>1
%     g = zeros(6,1);
% 
%   %W = Lt'\(Lt\eye(Var_target))-alpha_t*alpha_t';      % precompute for convenience
%   %Mt= -alpha_t;
%     W =Lt'\(Lt\eye(Var_target));
%     count = 1;
% 
%   for zz = 1 : 6
% %     if ( z <= nlf )                          % Gradient wrt  Kf
% %       [o p] = pos2ind_tri_inchol(z,m,irank); % determines row and column
% %       J = zeros(m,m); J(o,p) = 1;
% %       Val = J*Lf' + Lf*J';      
% %       dK = Val(ind_kf,ind_kf).*Kx11(ind_kx,ind_kx);
%     if zz ==1
%       Kn = L'\(L\eye(Var_source));
%       dK = -Kx12'*Kn*K12-K21*Kn*Kx12;
%       dMu = Kx12'*Kn*source_y;
%     end
%     if zz ==2
%         dK = Kx22;
%         dMu = zeros(Var_target,1);
%     end
%     
%     if zz >= 3 && zz<=4           % Gradient wrt parameters of Kx
%      
%       dKtt = adptivecovSEard(theta_x, target_x, zz)*theta_lamda2;
% 
%       dKst = adptivecovSEard(theta_x, x_source, zz,target_x)*theta_lamda1;
% 
%       Kn = L'\(L\eye(Var_source));
% 
%       dKss = adptivecovSEard(theta_x, x_source, zz);
% 
%       dKn = -Kn*dKss*Kn;
% 
%       dK =dKtt-dKst'*Kn*K12-K21*dKn*K12-K21*Kn*dKst;
%       
%       dMu =dKst'*Kn*source_y+K21*dKn*source_y;
%     end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       if ( z == (ltheta_x+1) )         % Gradient wrt theta_e variances
% %       
% %        Kn = L'\(L\eye(Var_source));
% %        dE = 2*exp(-theta_e)/((1+exp(-theta_e))^2);
% %        dK = Kx12'*dE*Kn*K12+K21*Kn*dE*Kx12;
% %       dMu = -Kx12'*dE*Kn*source_y;
% %      end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %       
% %     if ( zz == (ltheta_x+1) )         % Gradient wrt theta_b variances
% %      
% %       Kn = L'\(L\eye(Var_source));
% %       dK = Kx12'*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kn*K12+K21*Kn*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kx12;
% %       dMu = -Kx12'*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kn*source_y;
% %     end
% %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %   
% %     if ( zz == (ltheta_x+2) )         % Gradient wrt theta_mu variances
% %       Kn = L'\(L\eye(Var_source));
% %       dK = Kx12'*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kn*K12-K21*Kn*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kx12;
% %       dMu = -Kx12'*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kn*source_y;
% % 
% %     end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       
%     if ( zz == 5 )         % Gradient wrt Noise source variances
%       Kn = L'\(L\eye(Var_source));
%       
%       dKn = -Kn*2*noise_source*eye(Var_source)*Kn;
%       dK = -K21*dKn*K12;
%       dMu = K21*dKn*source_y;
% 
%     end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       
%     if ( zz == 6 )         % Gradient wrt Noise target variances
%  
%       dK = 2*noise_target*eye(Var_target);  
%       dMu = zeros(Var_target,1);
%       
% 
%      end % endif z
%       %W = Lt'\(Lt\eye(Var_target))-alpha_t*alpha_t';      % precompute for convenience
%   %Mt= -alpha_t;
%   
%  %gradnl(count) =  sum(sum(W.*dK,2),1)/2-dMu'*alpha_t/2;
%     g(count)=(trace(W*dK)-alpha_t'*dK*alpha_t-dMu'*alpha_t)/2;
%   
%     count = count + 1;
% 
%   end % end for derivarives   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else 
%     f = Inf;
% end
end
