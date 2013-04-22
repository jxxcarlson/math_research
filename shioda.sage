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


  One can list the characters for the (p,p) cohomology, along
  with a flag that signals whether chacacater space is spanned
  by a Hodge cycle:

    print(hodge_chars(2,5))==> 
      [[(3, 3, 2, 2) mod 5, True],
      [(3, 3, 3, 1) mod 5, False],
      [(4, 2, 2, 2) mod 5, False],
      [(4, 3, 2, 1) mod 5, True],
      [(4, 4, 1, 1) mod 5, True]]
  """

  def __init__(self, vector, modulus): # self, tuple, integer
    self.vector = vector
    self.modulus = modulus

  def __repr__(self):
    return "%s mod %d" % (self.vector, self.modulus)

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__dict__ == other.__dict__
    else:
      return False

  def __ne__(self, other):
    return not self.__eq__(other)

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

  def multiplicity(self, type):
    m = list(self.vector)
    if type == "s": # under action of symmetric group
      return Permutations(m).cardinality()
    if type == "g": # under action of galois group
      return self.galois_orbit().length()

  def hodge(self):
    """Return True if self is the character of a Hodge class."""
    w, d = self.vector, self.modulus
    a = self.q1()
    hodge_type = self.ht()
    if hodge_type[0] != hodge_type[1]:
      return "-"
    value = True
    for k in range(1,d-1):
      if gcd(k,d) == 1:
        c = self.mul(k)
        if c.q1() != a:
          value = False
    return value

  def galois_orbit(self):
    w = self.dim()
    value = []
    for k in range(1, self.modulus):
      if gcd(k, self.modulus) == 1:
        char = self.mul(k)
        value.append(char)
    value = unique(value)
    return USet(value)

class USet:
  
  def __init__(self, elements):
    self.elements = elements

  def __iter__(self):
    return iter(self.elements)

  def __repr__(self):
    value = "\n"
    for element in self:
      value = value + "%s\n" % element
    return value

  def length(self):
    return len(self.elements)

  def __eq__(self, other):

    if self.length() != other.length():
      return False

    value = True

    for element1 in self.elements:

      inner_value = False

      for element2 in other.elements:
        if element1 == element2:
          inner_value = True

      if inner_value == False:
        return False

    return value


def unique(L):
  '''Return sublist consisting of unique elements of L.'''
  M = []
  for item in L:
    if not item in M:
      M.append(item)
  return M

def in_related(item, L, relation):
  '''Return True if item is related to an item of L.'''
  for item2 in L:
    if relation(item, item2):
      return True
  return False

def unique2(L, relation):
  '''Return sublist of one element of L for each equivalence class of elements in L.'''
  M = []
  for item in L:
    if not in_related(item, M, relation):
      M.append(item)
  return M

class Character_augmented:
  
  def __init__(self, char, multiplicity, hodge): # self, Character, integer, integer
    self.char = char
    self.multiplicity = multiplicity
    self.hodge = hodge

  def __repr__(self):
    if self.hodge == -1:
      return "%s, %3d" % (self.char, self.multiplicity)
    if self.hodge == False:
      return "%s, %3d, False" % (self.char, self.multiplicity)
    elif self.hodge == True:
      return "%s, %3d, True" % (self.char, self.multiplicity)
    else:
      return "ERROR"

  def galois_orbit():
    return self.char.galois_orbit()      

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
  """Return list of monomials of degree d in n+2 variable of given weight."""
  degrees = [d for k in range(0,n+2)]
  value = []
  for monomial in mrange(degrees):
    if sum(monomial) == weight:
      value.append(monomial)
  return value

def monomial_classes(n, d, weight):
  """Return list of equivalence classes under the action of the symmetric 
  group of monomials of degree d in n+2 variable of given weight."""
  monoms = monomials(n, d, weight)
  value = []
  for m in monoms:
    if is_ordered(m):
      value.append(m)
  return value

def characters(n, d, weight):
  """Return list of chacters for the Fermat varietys of degree d 
  in n+2 variable of given weight."""
  monoms = monomials(n, d, weight)
  monoms = filter(has_positive_degrees, monoms)
  chars = []
  for m in monoms:
    chars.append(Character(tuple(m), d))
  return chars

def character_classes(n, d, weight):
  """Return list of equivalence classes under the action of the symmetric 
  group of monomials of degree d in n+2 variable of given weight."""
  monoms = monomial_classes(n, d, weight)
  monoms = filter(has_positive_degrees, monoms)
  chars = []
  for m in monoms:
    char = Character(tuple(m), d)
    """
    multiplicity = Permutations(m).cardinality()
    n = char.dim()
    p = n//2
    q = n - p
    if d*(q+1) == weight:
      is_hodge = char.hodge()
    else:
      is_hodge = -1
    ca = Character_augmented(char, Permutations(m).cardinality(), is_hodge)
    chars.append(ca)
    """
    chars.append(char)
  return chars

def hodge_chars(n,d, p):
  """Return list of classes characters for H^{p,q}"""
  q = n - p
  weight = (q+1)*d
  chars = character_classes(n,d,weight)
  return chars

def galois_orbits(n,d,p):
  """Return list of orbits given list of augmented character classes."""
  go = lambda x: x.galois_orbit()
  orbits = map(go, hodge_chars(n,d,p))
  orbits = unique(orbits)
  orbits = unique2(orbits, go_equivalent)
  return orbits



def go(n,d):
  '''Return one s-g equivalence class of character orbits for the Hodge structure of X(n,d).'''
  oo = []
  for q in range(0,n//2+1):
    p = n - q
    oo = oo + galois_orbits(n,d,p)
  oo = unique2(oo, go_equivalent)

  adjoin_info = lambda x: (galois_orbit_multiplicity(x), adjoin_hodge_types(x))
  oo = map(adjoin_info, oo)
  return oo

def print_orbit(o):
  print("")
  multiplicity = o[0]
  orbit = o[1]
  print multiplicity
  for ch in orbit:
    print(ch)

def pgo(n,d):
  '''Print orbit types, multiplicities, Hodge types for hypersurface
  of dimension n and degree d.'''
  oo = go(n,d)
  dim = 0
  index = 0
  for o in oo:
    print("")
    index += 1
    multiplicity = o[0]
    orbit = o[1]
    dim_orbit = len(orbit)
    dim += multiplicity*dim_orbit
    print("%d: dim = %d, mult = %d, total dim = %d" % (index, dim_orbit,  multiplicity, multiplicity*dim_orbit))
    for ch in orbit:
      print("%s %s" % (ch, ch[0].hodge()))
  print("")
  print("dimension = %d" % dim)
  print("")

def adjoin_hodge_types(o): ## x is a galois orbit
  o = list(o)
  oo = []
  for char in o:
    oo.append((char, char.ht()))
  return oo

def permutation_equivalent(perm1, perm2):
  if len(perm1) != len(perm2):
    return False
  P = Permutations(perm2)
  for p in P:
    if perm1 == p:
      return True
  return False

def permutation_in_list(p, L):
  """L is a liset of permutations.  If p is equivalent to an element of L, return True."""
  for q in L:
    if permutation_equivalent(p, q):
      return True
  return False

def count_equivalent_permutations(p, L):
  if type(p) == tuple:
    p = list(p)
  if type(L[0]) == tuple:
    L = map(lambda x: tuple(x), L)
  count = 0
  for q in L:
    if permutation_equivalent(p,q):
      count = count + 1
  return count

def ch2monom(ch):
  return ch.vector

def go_equivalent(o1, o2):
  o1 = list(o1)
  o2 = list(o2)
  m1 = list(o1[0].vector)
  M2 = map(lambda x: list(x.vector), o2)
  return permutation_in_list(m1, M2)

def galois_orbit_multiplicity(orbit):
  orbit = list(orbit)
  element = orbit[0]
  p = element.vector
  pp = map(lambda x: x.vector, orbit)
  inner_multiplicity = count_equivalent_permutations(p,pp)
  return element.multiplicity("s")/inner_multiplicity


def print_galois_orbits(n,d,p):
  orbits = galois_orbits(n,d,p)
  dim = 0
  for orbit in orbits:
    orb, mul = orbit
    dim += mul*len(orb)
    print(mul)
    print_orbit(orb)
    print("")
  print("Sum of dimensions = %d\n" % dim)



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

  sage: print_list(character_classes(2,5,10))
   ((3, 3, 2, 2) mod 5, 6)
   ((3, 3, 3, 1) mod 5, 4)
   ((4, 2, 2, 2) mod 5, 4)
   ((4, 3, 2, 1) mod 5, 24)
   ((4, 4, 1, 1) mod 5, 6)

   sage: print_list(hodge_class)chars(2,5))
   ((3, 3, 2, 2) mod 5, 6, True)
   ((3, 3, 3, 1) mod 5, 4, False)
   ((4, 2, 2, 2) mod 5, 4, False)
   ((4, 3, 2, 1) mod 5, 24, True)
   ((4, 4, 1, 1) mod 5, 6, True)
  """