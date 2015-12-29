function y = h_2_exponential(lx,lw)

% h_2_exponential(lx,lw)
% Calculates the entropy of y = x + w, where x and w are 
% two independent, arbitrarily-distributed exponential random variables.
%
% lx and lw are the (lambda) parameters of x and w, respectively
%
% lx and lw can be either scalars or matrices, but they must have
% the same dimensions

eg = double(eulergamma);

% ensure that lx < lw everywhere (makes no difference to result)
foo = sort([lx;lw],1);
lx = foo(1,:);
lw = foo(2,:);

% entropy is calculated differently depending on whether lx == lw
mask = (lx == lw);

% equal part
yeq = mask.*(2 - psi(2) - log(lw));

% not equal part
% force not equal to prevent numerical problems
lw = lw + mask;
yneq = (1+eg+psi(1+lx./(lw-lx))-log(lw.*lx./(lw-lx)))./((lw-lx)./lw);
yneq = yneq - ((lx./lw)+eg+psi(1+lw./(lw-lx))-log(lw.*lx./(lw-lx)))./((lw-lx)./lx);
yneq = (1-mask).*yneq;

y = yeq + yneq;