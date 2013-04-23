# An experiment to see how many totall real quintic fields of discriminant < N
# satisfy Nori's crtierion

# return the class number of the field defined by F:
def classNumber(F):
  R.<x> = PolynomialRing(QQ)
  f = R(F)
  K.<a> = NumberField(f)
  return K.class_number()

# Test to the field element x is epsilon-positive:

def  epsilon_positive(x):
  pos_minus_neg = sum( sgn(tx) for tx in x.complex_embeddings() )
  return ( pos_minus_neg == 3 )

    
def test0(F):
  R.<x> = PolynomialRing(QQ)
  f = R(F)
  K0.<a> = NumberField(f);
  D0 = K0.different(); d0 = D0.gens_reduced()[0]
  u = K0.units()

  i, j, k, l = 0, 0, 0, 0
  count = 0
  for i in range(0,2):
    for j in range(0,2):
      for k in range(0,2):
        for l in range(0,2):
          for m in range(0,2):
            sign = (-1)^i
            dd = sign*d0*u[0]^j*u[1]^k*u[2]^l*u[3]^m
            if epsilon_positive(dd) == True:
              # print "  ", dd.complex_embeddings()
              return True
  return False

def test(F):

  R.<x> = PolynomialRing(QQ)
  f = R(F)
  K0.<a> = NumberField(f);
  D0 = K0.different(); d0 = D0.gens_reduced()[0]
  u = K0.units()

  exponents = [2,2,2,2,2]
  count = 0
  for e in mrange(exponents):
    i,j,k,l,m = e
    sign = (-1)^i
    dd = sign*d0*u[0]^j*u[1]^k*u[2]^l*u[3]^m
    if epsilon_positive(dd) == True:
      # print "  ", dd.complex_embeddings()
      return True
  return False


# Apply the test to each field of discriminant < N

def testFields(N):
  print "Enumerating totally real fields of discriminant <", N
  print "This may take a while."

  TRF = enumerate_totallyreal_fields_all(5, N)
  print "Number of totally real fields:", len(TRF)

  n = 1
  nNonUFD = 0
  nPass = 0
  ratio = 0;

  for field in TRF:
    discr, G = field
    cn = classNumber(G)
    result = False
    if (cn == 1):
      result = test(G)
      if (result):
        nPass = nPass + 1
    else:
      nNonUFD = nNonUFD + 1
    ratio = nPass/(1.0*n)
    print "%5d %4d %s" % (n, discr, result)
    n = n + 1

  print "Summary:"
  print "  Number of fields:", len(TRF)
  print "  Number of fields of class number > 1:", nNonUFD
  print "  Number of fields which satisfy the criterion:", nPass
  print "  Ratio:", ratio
  
if __name__ == "__main__":
  testFields(10^6)