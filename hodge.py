"""
hodge.py: computation of Hodge numbers, etc., for 
hypersurfaces in projective space.
"""

def foo(n,k):
  """Return product n(n-1)...(k+1)k.
  >>> foo(4,1)
  24
  >>> foo(4,3)
  12
  """
  m = 1
  for i in range(k, n+1):
    m = m*i
  return m


def binomial(n,k):
  """Return binomial coefficient (n,k).
  >>> binomial(4,2)
  6
  """
  return foo(n,n-k+1)//foo(k,1)



def dim(d,n):
  """dimension of the space of homogeneous forms 
  of degree d in n+1 variables.
  >>> dim(3,2)
  10
  """
  return binomial(n+d, n)


def dim_moduli(d,n):
  """Dimension of space of moduli of hypersurfaces
  of degree d and dimension n.
  >>> dim_moduli(3,1)
  1
  >>> dim_moduli(3,2)
  4
  """
  return dim(d,n+1) - (n+2)*(n+2)


  #####################################################################   

   # j(d,n,r) = dimension of degree r component of the Jacobian ring       
   #            for a homogeneous form of degree d in n +1 variables.      
   #            The form is assumed to define a smooth variety.            

def dim_jacobian_ring(d,n,r):
  """Dimension of jacobian ring in degree r for a polynomial
  of degree d in n+1 variables.
  >>> dim_jacobian_ring(3,2,3)
  1
  >>> d, n = 5, 5
  >>> socle = (n+1)*(d-2)
  >>> dim_jacobian_ring(d, n, socle)
  1
  >>> k = socle//40
  >>> dim_jacobian_ring(d, n, k) == dim_jacobian_ring(d, n, socle - k)
  True
  >>> dim_jacobian_ring(d, n, d) == dim_moduli(d,n-1)
  True
  """
  dimj = 0
  i = 0
  rr = r
  while rr >= 0:
    term = (-1)**i * binomial(n+1, i) * dim(rr, n)
    dimj += term
    # print(i, rr, term, dimj)
    i = i + 1
    rr = rr - (d-1)
  return dimj

j = dim_jacobian_ring


def h(p,q,d,i=-1):
  """h(p,q,d) = dimension of H^{p,q} of hypersurface of degree d
   and dimension p+q.

   h(p,q,d,i) = same, but for i-th eigenspace  in case of cyclic
   hypersurface.

   >>> h(1,0,3)
   1
   >>> h(2,0,4)
   1
   >>> h(1,1,4)
   19
   >>> h(2,0,5)
   4
   >>> h(1,1,5)
   44
   >>> h(2,0,5,1)
   3
   >>> h(1,1,5,1)
   10
   >>> h(1,1,5,2)
   12
   >>> h(1,1,5,3)
   12
   >>> h(1,1,5,4)
   10
   >>> h(2,1,3,1)
   4
   >>> h(2,1,3,2)
   1
   >>> h(1,2,3,1)
   1
   >>> h(1,2,3,2)
   4
   
   """ 
  if i > -1:
    return j(d, p+q, (q+1)*d - (p+q+1) - i)
  else:
    return j(d, p+q+1, (q+1)*d - (p+q+2))


def betti(d,n):
  """Dimension of middle cohomology group
  for a hypersurface of degree d and dimenstion n.
  >>> betti(3,1)
  2
  >>> betti(4,1)
  6
  >>> betti(3,2)
  7
  >>> betti(4,2)
  22
  """
  b = 0
  for k in range(0,n+1):
    b += h(k,n-k,d)
  if n % 2 == 1:
    return b
  else:
    return b + 1


def euler_char(d,n):
  """Euler characteristic of a hypersuerface of degree n
  and dimension n.
  >>> euler_char(3,1)
  0
  >>> euler_char(3,2)
  9
  """
  if n == 0:
    return d
  else:
    return d*(n+1 - euler_char(d,n-1)) + euler_char(d,n-1)

ech = euler_char

def betti2(d):
  """Second Betti numbrer of an algebraic surface of 
  degree d in P^3 via different method than betti(d,3).
  >>> betti2(4) == betti(4,2)
  True
  >>> betti2(5) == betti(5,2)
  True
  """
  return ech(d,2) - 2

def hh(p,q,d):
  """Alternative method for computing h(p,q,d) ---
  sum up the h(p,q,d,i}.
  >>> p,q,d = 1,1,5
  >>> h(p,q,d) == hh(p,q,d)
  True
  """
  hpq = 0
  for i in range(1,d):
    hpq += h(p,q,d,i)
  return hpq

def dim_so(n):
  """Dimwsion of special orthogonal group.
  >>> dim_so(3)
  3
  """
  return n*(n-1)//2

def dim_u(n):
  """Diimension of the unitary group.
  >>> dim_u(2)
  4
  """
  return n*n

def dim_perdom2(a,b,c):
  """Complex dimension of a period domain of weight two
  with Hodge numbers a,b,c
  >>> dim_perdom2(1,19,1)
  19
  """
  dim = dim_so(a+b+c) - dim_u(a) - dim_so(b)
  return dim//2

def dim_udom2(p,q):
  """Complex dimension of U(p,q)/U(p)xU(q)
  >>> dim_udom2(1,10)
  10
  """
  dim = dim_u(p+q) - dim_u(p) - dim_u(q)
  return dim//2

if __name__ == "__main__": 

  # SELF TEST (runs docstring at head of module):
  from doctest import testmod
  import hodge
  testmod(hodge)
