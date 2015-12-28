function y = h_2_exponential(lx,lw)

foo = sort([lx;lw],1);
lx = foo(1,:);
lw = foo(2,:);

if (lx == lw)
    y = 2 - psi(2) - log(lw);
else
    eg = double(eulergamma);
    y = (1+eg+psi(1+lx./(lw-lx))-log(lw.*lx./(lw-lx)))./((lw-lx)./lw);
    y = y - ((lx./lw)+eg+psi(1+lw./(lw-lx))-log(lw.*lx./(lw-lx)))./((lw-lx)./lx);
end
