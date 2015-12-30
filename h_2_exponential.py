import scipy as sp
import scipy.special as sps

def h_2_exponential(lx,lw):
  # h_2_exponential(lx,lw)
  # 
  # Calculates the entropy of y = x + w, where x and w are 
  # two independent, arbitrarily-distributed exponential random variables.
  #
  # Requires scipy.
  #
  # lx and lw are the (lambda) parameters of x and w, respectively
  #
  # lx and lw can be either scalars or scipy arrays, but they must have
  # the same dimensions
  foo = sp.sort(sp.vstack((lx,lw)),0)
  hy = [h_noteq(r[0],r[1]) if (r[0]!=r[1]) else h_eq(r[0]) for r in sp.transpose(foo)]
  return sp.array(hy)

def h_eq(lw):
  # calculates entropy when the parameters are equal
  return (2 - sp.special.psi(2) - sp.log(lw));

def h_noteq(lx,lw):
  # calculates entropy when the parameters are nonequal
  eg = sp.euler_gamma
  hy = (1+eg+sps.psi(1+lx/(lw-lx))-sp.log(lw*lx/(lw-lx)))/((lw-lx)/lw)
  hy = hy - ((lx/lw)+eg+sps.psi(1+lw/(lw-lx))-sp.log(lw*lx/(lw-lx)))/((lw-lx)/lx)
  return hy

