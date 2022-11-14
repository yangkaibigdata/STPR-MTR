function [nl, gradnl] = nmargl_mtgp(logtheta, logtheta_all, ...
                          covfunc_x, x_source, f_source, x_target_train, f_target_train, D, deriv_range)
% Marginal likelihood and its gradients for multi-task Gaussian Processes
% 
% [nl, gradnl] = nmargl_mtgp(logtheta, logtheta_all, covfunc_x, x, y,...
%                	      m, irank, nx, ind_kf, ind_kx, deriv_range)
%
% To be used in conjunction with Carl Rasmussen's minimize function
% and the gpml package http://www.gaussianprocess.org/gpml/code/matlab/doc/
%
% nl = nmargl_mtgp(logtheta, ...) Returns the marginal negative log-likelihood
% [nl gradnl] =  nmargl_mtgp(logtheta, ...) Returns the gradients wrt logtheta
%
% logtheta    : Column vector of initial values of parameters to be optimized
% logtheta_all: Vector of all parameters: [theta_lf; theta_x; sigma_l]
%                - theta_lf: the parameter vector of the
%                   cholesky decomposition of k_f
%                - theta_x: the parameters of K^x
%                - sigma_l: The log of the noise std deviations for each task
% covfunc_x   : Name of covariance function on input space x
% x           : Unique input points on all tasks 
% y           : Vector of target values
% m           : The number of tasks
% irank       : The rank of K^f 
% nx          : number of times each element of y has been observed 
%                usually nx(i)=1 unless the corresponding y is an average
% ind_kx      : Vector containing the indexes of the data-points in x
%                which each observation y corresponds to
% ind_kf      : Vector containing the indexes of the task to which
%                each observation y corresponds to
% deriv_range : The indices of the parameters in logtheta_all
%                to which each element in logtheta corresponds to
%

% Author: Edwin V. Bonilla
% Last update: 23/01/2011


% *** General settings here ****
config = get_mtgp_config();
MIN_NOISE = config.MIN_NOISE;

% ******************************

if ischar(covfunc_x), covfunc_x = cellstr(covfunc_x); end % convert to cell if needed

%D = size(x,2); %  Dimensionality to be used by call to covfunc_x
n = length(f_target_train); 
m = length(f_source);
logtheta_all(deriv_range) = logtheta;
ltheta_x = D+1;     % number of parameters for input covariance

% nlf = irank*(2*m - irank +1)/2;        % number of parameters for Lf
% vlf = logtheta_all(1:nlf);             % parameters for Lf
% 
% theta_lf = vlf; 
% Lf = vec2lowtri_inchol(theta_lf,m,irank);
logtheta_all(ltheta_x+1)
logtheta_all(ltheta_x+2)
theta_x = logtheta_all(1:ltheta_x);                         % cov_x parameters
theta_b = abs(logtheta_all(ltheta_x+1));
theta_mu = abs(logtheta_all(ltheta_x+2));
%theta_e = logtheta_all(ltheta_x+1);
%sigma_source = exp(2*logtheta_all(ltheta_x+2));              % Noise parameters
%sigma_target = exp(2*logtheta_all(ltheta_x+3));                                        % Noise Matrix
sigma_source = exp(2*logtheta_all(ltheta_x+3));              % Noise parameters
sigma_target = exp(2*logtheta_all(ltheta_x+4));                                        % Noise Matrix
Var_source = size(x_source,1);
Var_target = size(x_target_train,1);
source_x = x_source;
target_x = x_target_train;
source_y = f_source;
target_y = f_target_train;

Kx11 = adptivecovSEard(theta_x, source_x);
K11 = Kx11+sigma_source*eye(Var_source);
[Kx22, Kx12] = adptivecovSEard( theta_x, source_x, target_x);
K22 = Kx22+sigma_target*eye(Var_target);
lamda = 2*(1/(1+theta_mu))^theta_b-1;
K12 = Kx12*lamda;
%K12 = Kx12*(2*(1/(1+theta_mu))^theta_b-1);
K21 = K12';


L =chol(K11)';
alpha = solve_chol(L',source_y);
%mu_target = sum(target_y)/Var_target;
mu_t = K21*alpha;
v   = L\K21';

C_t = K22-v'*v;

Lt = chol(C_t)';
y = target_y-mu_t;
alpha_t = solve_chol(Lt',y);

n=Var_source+Var_target;

% negative log-likelihood
nl = 0.5*y'*alpha_t + sum(log(diag(Lt))) + 0.5*n*log(2*pi);
%nl = -0.5*y'*alpha_t - sum(log(diag(Lt))) - 0.5*n*log(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 2)                      % If requested, its partial derivatives
  gradnl = zeros(size(logtheta));      % set the size of the derivative vector

  %W = Lt'\(Lt\eye(Var_target))-alpha_t*alpha_t';      % precompute for convenience
  %Mt= -alpha_t;
  W =Lt'\(Lt\eye(Var_target));
  count = 1;

  for zz = 1 : length(deriv_range) 
     z = deriv_range(zz);

%     if ( z <= nlf )                          % Gradient wrt  Kf
%       [o p] = pos2ind_tri_inchol(z,m,irank); % determines row and column
%       J = zeros(m,m); J(o,p) = 1;
%       Val = J*Lf' + Lf*J';      
%       dK = Val(ind_kf,ind_kf).*Kx11(ind_kx,ind_kx);

    if ( z <= ltheta_x )           % Gradient wrt parameters of Kx
     
      dKtt = adptivecovSEard(theta_x, target_x, z);

      dKst = adptivecovSEard(theta_x, x_source, z,target_x)*lamda;

      Kn = L'\(L\eye(Var_source));

      dKss = adptivecovSEard(theta_x, x_source, z);

      dKn = -Kn*dKss*Kn;

      dK =dKtt-dKst'*Kn*K12-K21*dKn*K12-K21*Kn*dKst;
      
      dMu =dKst'*alpha+K21*dKn*source_y;
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       if ( z == (ltheta_x+1) )         % Gradient wrt theta_e variances
%       
%        Kn = L'\(L\eye(Var_source));
%        dE = 2*exp(-theta_e)/((1+exp(-theta_e))^2);
%        dK = Kx12'*dE*Kn*K12+K21*Kn*dE*Kx12;
%       dMu = -Kx12'*dE*Kn*source_y;
%      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
    if ( z == (ltheta_x+1) )         % Gradient wrt theta_b variances
     
      Kn = L'\(L\eye(Var_source));
      dK = Kx12'*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kn*K12+K21*Kn*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kx12;
      dMu = -Kx12'*(2*log(1+theta_mu)*(1/(1+theta_mu))^theta_b)*Kn*source_y;
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
    if ( z == (ltheta_x+2) )         % Gradient wrt theta_mu variances
      Kn = L'\(L\eye(Var_source));
      dK = Kx12'*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kn*K12-K21*Kn*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kx12;
      dMu = -Kx12'*(2*theta_b*(1+theta_mu)^(-theta_b-1))*Kn*source_y;

    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
    if ( z == (ltheta_x+3) )         % Gradient wrt Noise source variances
      Kn = L'\(L\eye(Var_source));
      
      dKn = -Kn*2*sigma_source*eye(Var_source)*Kn;
      dK = -K21*dKn*K12;
      dMu = K21*dKn*source_y;

    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
    if ( z == (ltheta_x+4) )         % Gradient wrt Noise target variances
 
      dK = 2*sigma_target*eye(Var_target);  
      dMu = zeros(Var_target,1);
      

     end % endif z
      %W = Lt'\(Lt\eye(Var_target))-alpha_t*alpha_t';      % precompute for convenience
  %Mt= -alpha_t;
  
 %gradnl(count) =  sum(sum(W.*dK,2),1)/2-dMu'*alpha_t/2;
  gradnl(count)=(trace(W*dK)-alpha_t'*dK*alpha_t-dMu'*alpha_t)/2;
  
    count = count + 1;

  end % end for derivarives
  
end % end if nargout ==2





 

