"""
hodge_numbers.sage

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

  sage: moduli(6, [1,1,3]) = 3    -- double coversof P^1 branched at six points

  sage: moduli(6, [1, 2, 3]) = 1  -- elliptic curves in weighted projective space

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


# TESTS AND EXAMPLES:

r"""

sage: moduli(3, [1,1,1])     # moduli of cubic curves
  1

sage: moduli(3, 1)           # moduli of cubic curves
  1

sage: hodge(3,[1,1,1])       # hodge numbers of a cubic curve
  [1, 1]

sage: hodge(3,1)             # hodge numbers of a cubic curve
  [1, 1]

sage: hodge(4,[1,1,1,1])     # hodge numbers of a quartic surface
  [1, 19, 1]

sage: hodge(3,3,3)           # 3-sheeted cyclic cover of P^3
[0, 5, 5, 0]

sage: hodge(3,3,3,1)         # cyclic cubic threefold
  [0, 4, 1, 0]

sage: hodge(3,3,3,2)         # cyclic cubic threefold
  [0, 1, 4, 0]


"""

