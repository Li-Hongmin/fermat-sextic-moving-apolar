from itertools import product

# Exact second-order apolar transport for the existing raw t=-1 sum arc.
# A = K[s]/(s^3), K=Q(mu,a), Phi_12(mu)=0, a^3=2.
C.<mu> = CyclotomicField(12)
Za.<aa> = PolynomialRing(C)
K.<a> = C.extension(aa^3-2)
mu = K(mu)

TOP=(4,)*6

def comps(total,parts):
    if parts==1:
        yield (total,); return
    for first in range(total+1):
        for rest in comps(total-first,parts-1):
            yield (first,)+rest

def add(out,e,c):
    c=K(c);out[e]=out.get(e,K(0))+c
    if out[e]==0:del out[e]

def mul(A,B):
    out={}
    for e,c in A.items():
        for f,d in B.items():
            add(out,tuple(e[i]+f[i] for i in range(6)),c*d)
    return out

def base_aoki():
    out={};base=(1,0,0,2,0,0)
    qa=[((2,0,0,0,0,0),K(1)),((0,1,1,0,0,0),a)]
    cb=[((0,3,0,0,0,0),K(1)),((0,0,3,0,0,0),K(-1))]
    lc=[((0,0,0,0,j,4-j),mu^(4-j)) for j in range(5)]
    for (e,c),(f,d),(g,h) in product(qa,cb,lc):
        add(out,tuple(base[i]+e[i]+f[i]+g[i] for i in range(6)),c*d*h)
    return out

def permute_raw(A,sigma):
    out={}
    for e,c in A.items():
        f=[0]*6
        for old,new in enumerate(sigma):f[new]=e[old]
        add(out,tuple(f),c)
    return out

P0=base_aoki();P1=permute_raw(P0,(1,0,2,4,5,3));P=dict(P0)
for e,c in P1.items():add(P,e,-c) # geometric sum [S]+[S']

# h1 and h2 are exactly the coefficients selected by
# scratch/aoki_second_order_clean.sage.
h1={
 (0,0,2,0,4,0):-1,
 (0,0,2,1,3,0):-mu^2,
 (0,0,2,2,2,0):1-mu^2,
 (0,0,2,3,1,0):1,
 (0,0,2,4,0,0):mu^2,
 (1,1,0,0,4,0):a,
 (1,1,0,1,3,0):mu^2*a,
 (1,1,0,2,2,0):(mu^2-1)*a,
 (1,1,0,3,1,0):-a,
 (1,1,0,4,0,0):-mu^2*a,
}
h2={
 (0,0,4,0,0,2):K(1)/3*mu^2,
 (0,0,4,0,1,1):-K(1)/3*mu,
 (0,0,4,0,2,0):-K(1)/4,
 (0,0,4,1,0,1):-K(1)/3*mu^3,
 (0,0,4,2,0,0):(1-mu^2)/4,
 (0,3,1,0,0,2):K(7)/6*mu^2,
 (0,3,1,0,1,1):-K(2)/3*mu,
 (0,3,1,0,2,0):-K(1)/2,
 (1,1,2,0,0,2):-K(1)/2*mu^2*a,
 (1,1,2,1,0,1):-K(2)/3*mu^3*a,
 (1,1,2,2,0,0):(1-mu^2)/2*a,
}
h11=mul(h1,h1)

def reduced_pair(exp,pole):
    """Pair P with exact Fermat reduction from `pole` to pole 3."""
    steps=pole-3
    counts=[exp[i]//6 for i in range(6)]
    rem=[exp[i]-6*counts[i] for i in range(6)]
    if sum(counts)!=steps or sum(rem)!=12 or any(v==5 for v in rem):
        return K(0)
    denom=K(1)
    for current in range(pole,3,-1):denom*=K(6*(current-1))
    numer=K(1)
    for i in range(6):
        for j in range(counts[i]):numer*=K(exp[i]-5-6*j)
    comp=tuple(TOP[i]-rem[i] for i in range(6))
    if min(comp)<0:return K(0)
    return numer/denom*P.get(comp,K(0))

cache={}
def period_coeffs(A):
    """Coefficients of int_delta A Omega/F_s^3 through s^2."""
    if A in cache:return cache[A]
    z0=reduced_pair(A,3)
    z1=K(0)
    for e,c in h1.items():z1+=-3*c*reduced_pair(tuple(A[i]+e[i] for i in range(6)),4)
    z2=K(0)
    for e,c in h2.items():z2+=-3*c*reduced_pair(tuple(A[i]+e[i] for i in range(6)),4)
    for e,c in h11.items():z2+=6*c*reduced_pair(tuple(A[i]+e[i] for i in range(6)),5)
    cache[A]=(z0,z1,z2)
    return cache[A]

def _reduce_coefficient_row(row,basis):
    """Add one sparse row to an incremental echelon basis.

    Here ``row`` is a row of [M0|M1|M2].  Selecting a basis of these
    *coefficient triples* is legitimate over K[s]/(s^3): every discarded
    polynomial equation is the same constant K-linear combination of the
    retained polynomial equations at all three orders.
    """
    v=dict(row)
    while v:
        p=min(v)
        if p not in basis:
            lead=v[p]
            basis[p]={j:c/lead for j,c in v.items()}
            return True
        c=v[p]
        b=basis[p]
        for j,d in b.items():
            z=v.get(j,K(0))-c*d
            if z:
                v[j]=z
            elif j in v:
                del v[j]
    return False

def coefficient_basis_matrices(deg,predicate=None,label='full'):
    """Stream all monomial equations and retain only a joint row basis."""
    Bd=list(comps(deg,6));m=len(Bd)
    echelon={};selected=[];total=0;active=0
    for q in comps(12-deg,6):
        if predicate is not None and not predicate(q):
            continue
        total+=1
        rows=[{} for _ in range(3)]
        for j,g in enumerate(Bd):
            z=period_coeffs(tuple(q[k]+g[k] for k in range(6)))
            for r in range(3):
                if z[r]:rows[r][j]=z[r]
        packed={r*m+j:c for r in range(3) for j,c in rows[r].items()}
        if not packed:
            continue
        active+=1
        if _reduce_coefficient_row(packed,echelon):
            selected.append(rows)
    M=[matrix(K,len(selected),m,sparse=True) for _ in range(3)]
    for i,rows in enumerate(selected):
        for r in range(3):
            for j,c in rows[r].items():M[r][i,j]=c
    print('degree',deg,label,'streamed/active/joint-basis',total,active,len(selected),
          'matrix',len(selected),'x',m,'ranks',*[X.rank() for X in M],flush=True)
    return Bd,M

def derivative(A,i):
    out={}
    for e,c in A.items():
        if e[i]:
            f=list(e);f[i]-=1;add(out,tuple(f),c*e[i])
    return out

def shifted(poly,e):
    return {tuple(f[i]+e[i] for i in range(6)):c for f,c in poly.items()}

def apply_L(poly,order):
    return sum(c*period_coeffs(e)[order] for e,c in poly.items())

def block_invariants(M):
    M0,M1,M2=M;n=M0.nrows();m=M0.ncols()
    T2=block_matrix(K,[[M0,0],[M1,M0]],sparse=True)
    K2=T2.right_kernel_matrix()
    f2=K2.matrix_from_columns(range(m)).rank()
    T3=block_matrix(K,[[M0,0,0],[M1,M0,0],[M2,M1,M0]],sparse=True)
    K3=T3.right_kernel_matrix()
    f3=K3.matrix_from_columns(range(m)).rank()
    k0=M0.right_nullity()
    print('base_kernel/first_lift/second_lift',k0,f2,f3,flush=True)
    print('A-module = A^%d + (sA)^%d + (s^2A)^%d' % (f3,f2-f3,k0-f2),flush=True)
    return K2,K3,(k0,f2,f3)

def center_basis_crosscheck(deg,expected):
    _,Mr=coefficient_basis_matrices(deg,lambda e:max(e)<=4,'center')
    _,_,got=block_invariants(Mr)
    print('center Jacobian-basis invariant crosscheck',got,flush=True)
    assert got==expected

def poly(row,B,R):
    xs=R.gens()
    return sum(row[i]*prod(xs[j]^B[i][j] for j in range(6)) for i in range(len(B)))

R.<x0,x1,x2,x3,x4,x5>=PolynomialRing(K)

print('PERIOD TRANSPORT: L0(A)-3s L0(red4(Ah1))'
      '+s^2[-3L0(red4(Ah2))+6L0(red5(Ah1^2))]',flush=True)

# Strictness checks: the transported functional kills J(F_s), and the
# transported flat class stays of type (2,2), both modulo s^3.
Jfail=0;Jsample=None
dh1=[derivative(h1,i) for i in range(6)]
dh2=[derivative(h2,i) for i in range(6)]
for i in range(6):
    for r in comps(7,6):
        e0=list(r);e0[i]+=5;A0={tuple(e0):K(6)}
        A1=shifted(dh1[i],r);A2=shifted(dh2[i],r)
        vals=(apply_L(A0,0),
              apply_L(A0,1)+apply_L(A1,0),
              apply_L(A0,2)+apply_L(A1,1)+apply_L(A2,0))
        if any(vals):
            Jfail+=1
            if Jsample is None:Jsample=(i,r,vals)
print('Jacobian transport failures',Jfail,flush=True)
assert Jfail==0,Jsample

Hfail=0;Hsample=None
for q in comps(6,6):
    c1=sum(-2*c*reduced_pair(tuple(q[i]+e[i] for i in range(6)),3) for e,c in h1.items())
    c2=sum(-2*c*reduced_pair(tuple(q[i]+e[i] for i in range(6)),3) for e,c in h2.items())
    c2+=sum(3*c*reduced_pair(tuple(q[i]+e[i] for i in range(6)),4) for e,c in h11.items())
    if c1 or c2:
        Hfail+=1
        if Hsample is None:Hsample=(q,c1,c2)
print('Hodge-type transport failures',Hfail,flush=True)
assert Hfail==0,Hsample

# Fixed Fermat-basis coordinates of the transported period covector.
Pseries=[{} for _ in range(3)]
for e in product(range(5),repeat=6):
    if sum(e)!=12:continue
    comp=tuple(TOP[i]-e[i] for i in range(6))
    vals=period_coeffs(comp)
    for j,c in enumerate(vals):
        if c:Pseries[j][e]=c
assert Pseries[0]==P
print('period polynomial coefficient counts',*[len(z) for z in Pseries],flush=True)

B2,Mtwo=coefficient_basis_matrices(2)
K2two,K3two,inv2=block_invariants(Mtwo)
assert inv2==(1,1,1)
center_basis_crosscheck(2,inv2)
q0row=Mtwo[0].right_kernel_matrix().row(0)
q0=poly(q0row,B2,R)
print('central quadratic',q0,flush=True)

# Extract and normalize the unique free A-generator q(s).
m2=len(B2)
liftrow=next(row for row in K3two.rows() if any(row[j] for j in range(m2)))
g0=vector(K,liftrow[:m2]);g1=vector(K,liftrow[m2:2*m2]);g2=vector(K,liftrow[2*m2:])
idx=next(j for j in range(m2) if q0row[j])
scale=q0row[idx]/g0[idx]
g0*=scale;g1*=scale;g2*=scale
assert g0==q0row
assert Mtwo[0]*g0==0
assert Mtwo[1]*g0+Mtwo[0]*g1==0
assert Mtwo[2]*g0+Mtwo[1]*g1+Mtwo[0]*g2==0
print('lifted quadratic q1',poly(g1,B2,R),flush=True)
print('lifted quadratic q2',poly(g2,B2,R),flush=True)

# Test q(s)=l0(s)l1(s) through s^2.  In the normalization printed above,
# q0=(x4-mu*x5)(x3-mu^(-1)*x5).  First-order factorability means q1 lies
# in l0*S1+l1*S1.
l0=x4-mu*x5;l1=x3-mu^(-1)*x5
qbase=l0*l1
basevec=vector(K,[qbase.monomial_coefficient(prod(R.gen(j)^e[j] for j in range(6))) for e in B2])
facscale=K(1)
fg0=g0;fg1=g1;fg2=g2
assert fg0==basevec
cols=[]
for i in range(6):
    pp=R.gen(i)*l1
    cols.append(vector(K,[pp.monomial_coefficient(prod(R.gen(j)^e[j] for j in range(6))) for e in B2]))
for i in range(6):
    pp=l0*R.gen(i)
    cols.append(vector(K,[pp.monomial_coefficient(prod(R.gen(j)^e[j] for j in range(6))) for e in B2]))
Tfac=matrix(K,cols).transpose()
aug1=Tfac.augment(matrix(K,m2,1,list(fg1)))
first_factorable=(aug1.rank()==Tfac.rank())
print('factor tangent rank/augmented',Tfac.rank(),aug1.rank(),'first_factorable',first_factorable,flush=True)
assert not first_factorable
if first_factorable:
    w=Tfac.solve_right(fg1)
    u=sum(w[i]*R.gen(i) for i in range(6));v=sum(w[6+i]*R.gen(i) for i in range(6))
    residual=fg2-vector(K,[(u*v).monomial_coefficient(prod(R.gen(j)^e[j] for j in range(6))) for e in B2])
    aug2=Tfac.augment(matrix(K,m2,1,list(residual)))
    print('one-gauge second factor rank/augmented',Tfac.rank(),aug2.rank(),flush=True)

B3,Mthree=coefficient_basis_matrices(3)
K2three,K3three,inv3=block_invariants(Mthree)
center_basis_crosscheck(3,inv3)

# Degree-three primitive quotient at the central fibre.
K03=Mthree[0].right_kernel_matrix()
q0multiples=[]
for variable in R.gens():
    pp=q0*variable
    q0multiples.append(vector(K,[pp.monomial_coefficient(prod(R.gen(j)^e[j] for j in range(6))) for e in B3]))
Mult=matrix(K,q0multiples)
print('central cubic kernel / quadratic multiples / minimal-generator quotient',K03.nrows(),Mult.rank(),K03.stack(Mult).rank()-Mult.rank(),flush=True)
assert inv3==(10,10,10) and Mult.rank()==6

# Verify the transported quadratic generator really gives a free rank-six
# submodule q(s)S_1 inside Ann_3.  A unit 6-minor at s=0 then splits this
# inclusion over the local Artin ring, so the quotient is free of rank 4.
def multiply_by_variables(v,Bfrom,Bto):
    pos={e:i for i,e in enumerate(Bto)}
    C=matrix(K,len(Bto),6,sparse=True)
    for j,e in enumerate(Bfrom):
        if not v[j]:continue
        for i in range(6):
            f=list(e);f[i]+=1
            C[pos[tuple(f)],i]+=v[j]
    return C

C0=multiply_by_variables(g0,B2,B3)
C1=multiply_by_variables(g1,B2,B3)
C2=multiply_by_variables(g2,B2,B3)
assert Mthree[0]*C0==0
assert Mthree[0]*C1+Mthree[1]*C0==0
assert Mthree[0]*C2+Mthree[1]*C1+Mthree[2]*C0==0
assert C0.rank()==6
print('q(s)S1 inclusion strictness failures 0; central rank',C0.rank(),flush=True)
print('cubic minimal-generator quotient = A^4; truncated generators beta_1,2=1 beta_1,3=4 are flat',flush=True)

# A minimal certificate: ranks of the coefficient-block maps.
print('CERT ann2 invariants',inv2,flush=True)
print('CERT ann3 invariants',inv3,flush=True)
