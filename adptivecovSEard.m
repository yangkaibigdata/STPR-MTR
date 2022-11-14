function [A, B] = adptivecovSEard(logtheta, x, z, y)

% Squared Exponential covariance function with Automatic Relevance Detemination
% (ARD) distance measure. The covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2)
%
% where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
% D is the dimension of the input space and sf2 is the signal variance. The
% hyperparameters are:
%
% logtheta = [ log(ell_1)
%              log(ell_2)
%               .
%              log(sqrt(sf2)) ]%              log(ell_D)

%
% For more help on design of covariance functions, try "help covFunctions".
%
% (C) Copyright 2006 by Carl Edward Rasmussen (2006-03-24)

%if nargin == 0, A = D+1; end          % report number of parameters
%persistent K;    
%persistent L;
[n, D] = size(x);
 

%ell = exp(logtheta(1:D));                         % characteristic length scale
%sf2 = exp(2*logtheta(D+1));                                   % signal variance
%ell = exp(logtheta(1:D));                         % characteristic length scale
%sf2 = exp(2*logtheta(D+1));                                   % signal variance
ell = logtheta(1:D);                         % characteristic length scale
sf2 = 1;                                   % signal variance
%sf = logtheta(D+1);
if nargin == 2 && nargout==1
  A = sf2*exp(-sq_dist(diag(1./ell)*x')/2);
 
  B =0;
end
if nargout == 2                          % compute test set covariances
  A = sf2*exp(-sq_dist(diag(1./ell)*z')/2);
  %A = sf2*ones(size(z,1),1);
  B = sf2*exp(-sq_dist(diag(1./ell)*x',diag(1./ell)*z')/2); 
end
if nargin==3 && nargout==1                                                % compute derivative matrix
    K = sf2*exp(-sq_dist(diag(1./ell)*x')/2);
  if z <= D                                                      % length scale parameters
    
    A = K.*((sq_dist(x(:,z)'/ell(z)))/ell(z));
  %  A = K.*(sq_dist(x(:,z)'/ell(z)));
    B =0;
  else
    A = 2*sf*exp(-sq_dist(diag(1./ell)*x')/2);
    %A = 2*K;
    B =0;
  end
end
  
if nargin==4 && nargout==1
     L=sf2*exp(-sq_dist(diag(1./ell)*x',diag(1./ell)*y')/2);
     if z<=D
         A = L.*((sq_dist(x(:,z)'/ell(z),y(:,z)'/ell(z)))/ell(z));
         %  A = L.*(sq_dist(x(:,z)'/ell(z),y(:,z)'/ell(z)));
          B =0;
     
    else                                                        % magnitude parameter
     % A = 2*L;
      A = 2*sf*exp(-sq_dist(diag(1./ell)*x',diag(1./ell)*y')/2);
      B =0;
     end
end


return

