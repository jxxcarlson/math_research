"""
hodge_numbers.sage

Author: J. Carlson
Date:   March 28, 2013
URL:    www.math.utah.edu/~carlson

Compute hodge numbers of hypersurfaces
in weighted projective spaces, and of cyclic covers.
For the latter, compute also the hodge numbers of the 
eigenspaces of the covering automorphisms.

For examples, see the tests at the end of this file.
To load this package:

#  sage: load('hodge_numbers.sage')

Then to ask for help/documentation:

#  sage: hodge?
#  sage: moduli?

The functions 'hodge' and 'moduli' are the 
main functions in this module.

To run the tests, do this:

  $ sage -t whodge.sage

"""

var('t')

def poincare_f(d,q):
  return simplify((1-t^d)/(1-t^q))

def poincare_function(L):
  """
  Poincare polynomial for a regular sequence where
  L = [[d1, w1], [d2, w2], .. ] is a list of 
  degrees and weight
  """
  p = 1
  for pair in L:
    a, b = pair
    p = p*poincare_f(a,b)
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
  return poincare_function(jacobian_weights(d,W))
   
def  moduli(d,w):
  """
  moduli(d,w) is the dimension of the moduli space of hypersurfaces  of degree d 
  w can be a list of weights, e.g., moduli(6, [1,2,3]), or it can be a number,
  e.g., moduli(3,1).  

  In the first case, moduli(d,w) is the number of moduli of a hypersurface 
  of degree d in a weighted projective spaces with weights w.
  
  In the second case, it is the number of moduli of a hypersurface of degree d
  and dimension n.

  Examples:

  moduli(4,2) = 19                -- the number of moduli of quartic surfaces

  moduli(4, [1,1,1,1]) = 19       -- as above

  moduli(6, [1,1,3]) = 3    -- double coversof P^1 branched at six points

  moduli(6, [1, 2, 3]) = 1  -- elliptic curves in weighted projective space

  """
  if (type(w) == sage.rings.integer.Integer): # case of X of dim n in P^(n+1)
    w = [1 for r in range(0,w+2)]
  P = pj(d,w)
  T = P.taylor(t,0,d+1)
  return T.coefficient(t,d)

def HODGE(d,W,s=0):
  """
  HODGE(d,W): return vector of Hodge numbers for hypersurfce
  of degree d and  weight vector W.

  """
  P = pj(d,W)                     #  Poincare polynomial
  n = len(W)                      # number of homogeneous variables
  vdeg = sum(W)
  tdeg = (n+1)*d - vdeg + 1
  T = P.taylor(t,0,tdeg)
  H = [ ]
  if s > 0: # KLUDGE / HACK
    n = n + 1
  for q in range(0,n-1):
    c = T.coefficient(t,(q+1)*d - vdeg - s)  # Hodge number h(p,q)
    H.append(c)
  return H

def hodge(d, w, k=0, i=0):
  """
  hodge(d,n) = hodge vector for hypersurface of degree d and dimension n
  
  hodge(d,W) = hodge vector for hypersurface of degree d and weight vector W

  hodge(d,n,k) = hodge vector for k-to-1 cyclic cover of hypersurface 
  of degree d and dimension n.

  hodge(d,n,k,i) = hodge vector for k-to-1 cyclic cover of P^n branched  along
  a hypersurface of degree d.

  Examples:

  hodge(4,2) = [1,19,1]           -- K3 surface

  hodge(4,[1,1,1,1]) = [1,19,1]   -- as above

  hodge(6, [1, 2, 3]) = [1,1]     -- elliptic curves in P(1,2,3)

  hodge(6,2,2) = [1,19,1]         -- double covers of P^2 branched on a sextic curve

  hodge(3,3,3) = [0, 5, 5, 0]     -- cyclic cubic threefold: branched cover of P^3

  hodge(3,3,3,1) = [0, 4, 1, 0]   -- hodge numbers of an eigenspace of previous threefold

  hodge(3,3,3,2) = [0, 1, 4, 0]   -- the other eigenspace

  """
  if (type(w) == sage.rings.integer.Integer) & (k==0): # case of X of dim n in P^(n+1)
    w = [1 for r in range(0,w+2)]
   
  if (type(w) == sage.rings.integer.Integer) & (k>0) & (i==0): # k-sheeted cyclic cover of P^n
    w = [1 for r in range(0,w+1)]
    w.append(d//k)

  if (type(w) == sage.rings.integer.Integer) & (k>0) & (i>0): # k-sheeted cyclic cover of P^n
    w = [1 for r in range(0,w+1)]
    i = i*(d//k)
  
  return HODGE(d,w,i)




def betti(d,n):
  """Dimension of middle cohomology group                 
  for a hypersurface of degree d and dimenstion n.
  """
  b = sum(hodge(d,n))    # primitive middel cohomology
  if n % 2 == 1:
    return b
  else:
    return b + 1


def euler_char(d,n):
  """Euler characteristic of a hypersuerface of degree n
  and dimension n.                                      
  """
  if n == 0:
    return d
  else:
    return d*(n+1 - euler_char(d,n-1)) + euler_char(d,n-1)

def dim_so(n):
  """
  Dimension of special orthogonal group.
  """
  return n*(n-1)//2

def dim_u(n):
  """
  Dimension of the unitary group.                                                               
  """
  return n*n

def dim_perdom(hv):
  """Complex dimension of a period domain with given Hodge vector. """
  n = sum(hv)
  if len(hv) == 3:
    p,q,r = hv
    dim = dim_so(n) - dim_u(p) - dim_so(q)
    return dim//2
  else:
    return -1

def dim_horizontal(hv):
  """Complex dimension of horiziontal tangent bundle."""
  n = sum(hv)
  if len(hv) == 3:
    p,q,r = hv
    return p*q
  else:
    return -1

def dim_uherm(p,q):
  """Complex dimension of U(p,q)/U(p)xU(q)"""
  dim = dim_u(p+q) - dim_u(p) - dim_u(q)
  return dim//2

################################################################
#
#            Abbreviations
#
################################################################

ech = euler_char
dim_h = dim_horizontal

################################################################
#
# TESTS AND EXAMPLES:
#
################################################################


r"""
################################################################
#
#            Moduli
#
################################################################


sage: moduli(3, [1,1,1])        # moduli of cubic curves
  1

sage: moduli(3, 1)              # moduli of cubic curves
  1

sage: moduli(4,2)               #  moduli of quartic surfaces
  19

sage: moduli(4, [1,1,1,1])      # as above
  19

sage: moduli(6, [1,1,3])        # double coversof P^1 branched at six points
  3

sage: moduli(6, [1, 2, 3])      # elliptic curves in P(1,2,3)
  1


################################################################
#
#           Hodge numbers
#
################################################################

sage: hodge(3,[1,1,1])           # hodge numbers of a cubic curve
  [1, 1]

sage: hodge(3,1)                 # hodge numbers of a cubic curve
  [1, 1]

sage: hodge(4,[1,1,1,1])         # hodge numbers of a quartic surface
  [1, 19, 1]

sage: hodge(4,2)                 #  as above
  [1, 19, 1]

sage: hodge(6, [1, 2, 3])        # elliptic curves in P(1,2,3)
  [1, 1]

sage: hodge(6,2,2)               # double covers of P^2 branched on a sextic curve
  [1, 19, 1]

sage: hodge(3,3,3)               # 3-sheeted cyclic cover of P^3
[0, 5, 5, 0]

sage: hodge(3,3,3,1)             # cyclic cubic threefold
  [0, 4, 1, 0]

sage: hodge(3,3,3,2)             # cyclic cubic threefold
  [0, 1, 4, 0]

################################################################
#
#            Topogy
#
################################################################

sage: betti(3,1)                # cubic curve
  2

sage: betti(4,2)                # quartic surface
  22

sage: euler_char(3,1)           # cubic curve
  0

sage: ech(3,1)                  # the same (^)
  0

sage: ech(4,2)                   # 24 (XX: ALREADY?)
  24

################################################################
#
#            Period domains
#
################################################################

sage: dim_so(3)
  3

sage: dim_so(2)
  1

sage: dim_u(2)
  4

sage: dim_u(1)
  1

sage: dim_perdom([1,19,1])
  19

sage: dim_perdom([2,3,2])
  7

sage: dim_uherm(1,4)
  4

################################################################
#
#            Helper functions
#
################################################################

sage: poincare_f(4,2)           # Poincare function
  (t^4 - 1)/(t^2 - 1)

sage: jacobian_weights(4,[1,1,1])
  [[3, 1], [3, 1], [3, 1]]

sage: jacobian_weights(6,[1,2,3])
  [[5, 1], [4, 2], [3, 3]]

sage: pj(3,[1,1,1])
  t^6/(t - 1)^3 - 3*t^4/(t - 1)^3 + 3*t^2/(t - 1)^3 - 1/(t - 1)^3

sage: pj(6,[1,2,3])
  t^9/((t - 1)*(t^2 - 1)) - t^5/((t - 1)*(t^2 - 1)) - t^4/((t - 1)*(t^2 - 1)) + 1/((t - 1)*(t^2 - 1))

################################################################

"""

