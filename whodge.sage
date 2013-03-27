
# var('t,d,q')

var('t')

def poincare_function(d,q):
  return simplify((1-t^d)/(1-t^q))

def wpoincare_function(L):
  """
  Poincare polynomial for a regular sequence where
  L = [[d1, w1], [d2, w2], .. ] is a list of 
  degrees and weight
  """
  p = 1
  for pair in L:
    a, b = pair
    p = p*poincare_function(a,b)
  return expand(p)

def jacobian_weights(d,W):
  """ 
  jacobian_weights(d,W ) is the list of degrees and weights for 
  the regular sequence obtained from a weighted homogeneous
  form of degree where W is the list of weights
  """
  L = []
  for weight in W:
    L.append( [d-weight, weight] )
  return L
    
def pj(d,W):
  """
  pj(d,W) = Poincare polynomial for Jacobian ring of poly of degre d
  """
  return wpoincare_function(jacobian_weights(d,W))
   
def  wmoduli(W,d):
  """
  wmoduli( W, d ) is the dimension of the moduli space of hypersurfaces
  of degree d with weights W
  >>> wmoduli([1,1,1],3) # cubic curve
  1
  """
  P = pj(d,W)
  T = P.taylor(t,0,d+1)
  return T.coefficient(t,d)

def hodge(d,W):
  """
  hodge(d,W): return vector of Hodge numbers for hypersurfce
  of degree d and  weight.

  hodge(d,n): return vector of Hodge numbersf for a hypersurface
  of degree and dimension n

  >>> hodge(3,[1,1,1])       # cubic curve
  [1, 1]

  >>> hodge(4,[1,1,1,1])     # quartic surface
  [1, 19, 1]
  """
  if type(W) == sage.rings.integer.Integer:
    W = [1 for k in range(0,W+2)]
  P = pj(d,W)                     #  Poincare polynomial
  n = len(W)                      # number of homogeneous variables
  vdeg = sum(W)
  tdeg = (n+1)*d - vdeg + 1
  T = P.taylor(t,0,tdeg)
  H = [ ]
  for q in range(0,n-1):
    c = T.coefficient(t,(q+1)*d - vdeg)  # Hodge number h(p,q)
    H.append(c)
  return H

def chodge(k,d,m):
  """
  Return vector of Hodge numbers for a k-fold cyclic cover
  of P^m branched along a hypersurface of degree d, where
  k divides d.
  >>> chodge(3,3,3) # 3-sheeted cyclic cover of P^3
  [0, 5, 5, 0]
  """
  W = [1 for r in range(0,m+1)]
  W.append(d//k)
  P = pj(d,W)                     #  Poincare polynomial
  n = len(W)                      # number of homogeneous variables
  vdeg = sum(W)
  tdeg = (n+1)*d - vdeg + 1
  T = P.taylor(t,0,tdeg)
  H = [ ]
  for q in range(0,m+1):
    c = T.coefficient(t,(q+1)*d - vdeg)  # Hodge number h(p,q)
    H.append(c)
  return H


def echodge(k,d,m,i):
  """
  Return vector of Hodge numbers for the i-th eigenspace 
  of a k-fold cyclic cover of P^m branched along 
  a hypersurface of degree d, where k divides d.
  >>> echodge(3,3,3,1)  # cyclic cubic threefold
  [0, 4, 1, 0]
  >>> echodge(3,3,3,2)  # cyclic cubic threefold
  [0, 1, 4, 0]
  """
  W = [1 for r in range(0,m+1)]
  P = pj(d,W)                     #  Poincare polynomial
  n = len(W)                      # number of homogeneous variables
  vdeg = sum(W)
  tdeg = (n+1)*d - vdeg + 1
  T = P.taylor(t,0,tdeg)
  H = [ ]
  for q in range(0,m+1):
    c = T.coefficient(t,(q+1)*d - vdeg - i*(d//k) )  # Hodge number h(p,q)
    H.append(c)
  return H

if __name__ == "__main__":
  # SELF TEST (runs docstring at head of module):
  from doctest import testmod
  import whodge
  testmod(whodge)
