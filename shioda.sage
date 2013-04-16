class Character:
  """Class of characters for the finite abelian group G considered
  in Shioda's paper "The Hodge Conjecture for Fermat Varieties, Math.
  Ann. 245, 175-184 (1979).

  Test this file with 

    $ sage -t shioda.sage

  Example: let W be the contradtion of the Eulder vector field with 
  the "complex volume form" dV = dX_0 \wedge ... \wedge dX_{n+1}.  Let
  a = residue(MW/F^{q+1}), where M is a monomial in the X_i, and
  where F is the sum the d-th powers of the variables.  If the powers
  of the X_i in M are strictly less than d-1, then the residue
  is nonzero in cohomology.  Such residues jointly give a basis for 
  the Hodge components of the H^{p,q}. 

  Let F = X_0^5 + ... + X_3^5.  Let V be the variety defined by F = 0.

  Let M = X_0.

  Let r = residue(MW/F).  This is a class in H^{2,0)(V).  It has
  character (2,1,1,,1), where the entries are viewed in Z/5.
  In the present program, we construct the character associated 
  with the residue r by

    c0 = Character((2,1,1,1),5)

  This vector (2, 1, 1, 1) should be thought of as the sum
  of (1, 0, 0, 0), which is the charcter for X_0, and
  (1, 1, 1, 1), which is the character for W.  Since F is
  homogeneous of degree 5, its charater vector is (0,0,0,0).

  The Hodge type for the character is computed by 

    c0.ht()

  Now consider the monomial X_0^3X_1^3, with character
  (3, 3, 0, 0).  The character for the associated residue is
  (4, 4, 1, 1).  We find:

    c1 = Character((4,4,1,1), 5) 
    c1.ht() ==> (1,1)
    c2.hodge() ==> True

  The last line says the the residue with character c1 is
  a Hodge class.

  Here is a complete list of characters for the cohomology
  of the Fermat quintic surface:

    H^{2,0}:  (2,1,1,1)
    H^{1,1}:  (4,4,1,1), (4,3,2,1), (4,2,2,2), (3,3,3,1), (3,3,1,1) 
    H^{0,2):  (4,4,4,3)

  where H^{1,1} refers to primitive (1,1) cohomology.  One finds
  that for H^{1,1}, that the first three weight spaces consist
  of Hodge classes, while there are now Hodge classes in the 
  last to weight spaces.

  """

  def __init__(self, vector, modulus): # self, tuple, integer
    self.vector = vector
    self.modulus = modulus

  def __repr__(self):
    return "%s mod %d" % (self.vector, self.modulus)


  def dim(self):
    """Dimension of the underlying hypersurface."""
    return len(self.vector) - 2;  


  def q1(self):
    """If self is character of a class of type (p,q), return q+1."""
    w, d = self.vector, self.modulus
    value = 0
    for k in range(0, len(w)):
      a = w[k] % d
      value += a/d
    return value

  def ht(self):
    q = self.q1()
    n = self.dim()
    q = q-1
    p = n - q
    return (p,q)


  def mul(self, scalar):
    """Return self (character) multiplied by a scalar."""
    w, d = self.vector, self.modulus
    value = ()
    for k in range(0, len(w)):
      x = scalar * w[k] % d
      value = value + (x,)
    return Character(value, d)

  def test(self):
    w, d = self.vector, self.modulus
    print("%4d: %d %s" % (0, self.q1(), w))
    for k in range(1,d-1):
      print("%4d: %d -- %s" % (k, self.q1(), self.mul(k)))

  def hodge(self):
    """Return True if self is the character of a Hodge class."""
    w, d = self.vector, self.modulus
    a = self.q1()
    value = True
    for k in range(1,d-1):
      c = self.mul(k)
      if c.q1() != a:
        value = False
      # print(c)
    return value

def is_ordered(vector):
  value = True
  for i in range(1, len(vector)):
    if vector[i] > vector[i-1]:
      value = False
  return value

def has_positive_degrees(vector):
  value = True
  for i in range(0, len(vector)):
    if vector[i] < 1:
      value = False
  return value

def monomials(n, d, weight):
  degrees = [d for k in range(0,n+2)]
  value = []
  for monomial in mrange(degrees):
    if sum(monomial) == weight:
      value.append(monomial)
  return value

def monomial_classes(n, d, weight):
  monoms = monomials(n, d, weight)
  value = []
  for m in monoms:
    if is_ordered(m):
      value.append(m)
  return value

def characters(n, d, weight):
  monoms = monomials(n, d, weight)
  monoms = filter(has_positive_degrees, monoms)
  chars = []
  for m in monoms:
    chars.append(Character(tuple(m), d))
  return chars

def character_classes(n, d, weight):
  monoms = monomial_classes(n, d, weight)
  monoms = filter(has_positive_degrees, monoms)
  chars = []
  for m in monoms:
    chars.append(Character(tuple(m), d))
  return chars

def hodge_chars(n,d):
  if n % 2 == 1:
    return []
  p = n//2
  weight = (p+1)*d
  chars = character_classes(n,d,weight)
  value = []
  for char in chars:
    value.append([char, char.hodge()])
  return value

def print_list(L):
  for item in L:
    print(item)




  r"""
  sage: c1 = Character((4,4,1,1), 5)
  sage: c1
    (4, 4, 1, 1) mod 5

  sage: c1.q1()
    2
  
  sage: c1.dim()
    2

  sage: c1.ht()
    (1, 1)

  sage: c1.mul(2)
    (3, 3, 2, 2) mod 5

  sage: c1.hodge()
    True

  sage: c0 = Character((2,1,1,1),5)
  sage: c0.q1()
    1

  sage: c0.hodge()
    False

  sage: c = Character((1,0,0,0), 5)
  sage: c.ht()
    (14/5, -4/5)

  sage: foo = hodge_chars(2,5)
  sage: foo
    [[(3, 3, 2, 2) mod 5, True],
    [(3, 3, 3, 1) mod 5, False],
    [(4, 2, 2, 2) mod 5, False],
    [(4, 3, 2, 1) mod 5, True],
    [(4, 4, 1, 1) mod 5, True]]
  """