from itertools import product

# Exact coefficient field K = Q(mu,a), mu^4-mu^2+1=0, a^3=2.
C.<mu> = CyclotomicField(12)
Za.<aa> = PolynomialRing(C)
K.<a> = C.extension(aa^3 - 2)
mu = K(mu)

CAP = 4
SOCLE = (4,)*6
B6 = [e for e in product(range(5), repeat=6) if sum(e) == 6]
I6 = {e:i for i,e in enumerate(B6)}

def addterm(out, e, c):
    out[e] = out.get(e,K(0)) + c
    if out[e] == 0:
        del out[e]

def base_aoki():
    out = {}
    base = (1,0,0,2,0,0)
    fa = [((2,0,0,0,0,0),K(1)), ((0,1,1,0,0,0),a)]
    fb = [((0,3,0,0,0,0),K(1)), ((0,0,3,0,0,0),K(-1))]
    fc = [((0,0,0,0,j,4-j),mu^(4-j)) for j in range(5)]
    for (ea,ca),(eb,cb),(ec,cc) in product(fa,fb,fc):
        e=tuple(base[i]+ea[i]+eb[i]+ec[i] for i in range(6))
        addterm(out,e,ca*cb*cc)
    return out

def permute_variables(poly, sigma):
    out={}
    for e,c in poly.items():
        ee=[0]*6
        for old,new in enumerate(sigma):
            ee[new]=e[old]
        addterm(out,tuple(ee),c)
    return out

def rho(poly):
    M=matrix(K,426,426,sparse=True)
    for i,left in enumerate(B6):
        for e,c in poly.items():
            right=tuple(4-left[j]-e[j] for j in range(6))
            if min(right)>=0 and max(right)<=4 and right in I6:
                M[i,I6[right]] += c
    return M

P=base_aoki()
Q=permute_variables(P,(1,0,2,4,5,3))
A=rho(P); B=rho(Q)
print('global ranks',A.rank(),B.rank(),A.stack(B).rank())

# Simultaneous sparsity components.
adj=[set() for _ in range(426)]
for M in (A,B):
    for (i,j),c in M.dict().items():
        adj[i].add(j);adj[j].add(i)
seen=set(); comps=[]
for i in range(426):
    if i in seen: continue
    seen.add(i); todo=[i]; comp=[]
    while todo:
        u=todo.pop();comp.append(u)
        for v in adj[u]:
            if v not in seen:
                seen.add(v);todo.append(v)
    comps.append(sorted(comp))

sing=None
for comp in comps:
    AA=A.matrix_from_rows_and_columns(comp,comp)
    BB=B.matrix_from_rows_and_columns(comp,comp)
    ranks=(len(comp),AA.rank(),BB.rank(),AA.stack(BB).rank(),(AA+BB).rank())
    if len(comp)>1: print('component',ranks)
    if ranks[:4] == (60,12,12,23): sing=(comp,AA,BB)
assert sing is not None
comp,AA,BB=sing

# Quotient by the simultaneous radical.  If N is an echelon basis of the
# common radical, standard vectors in its non-pivot coordinates form a
# complement.
N=AA.stack(BB).right_kernel_matrix()
assert N.nrows()==37
pivots=N.echelon_form().pivots()
free=[j for j in range(60) if j not in pivots]
assert len(free)==23
U=matrix(K,60,23,lambda i,j: K(1) if i==free[j] else K(0))
assert N.stack(U.transpose()).rank()==60
Aq=U.transpose()*AA*U
Bq=U.transpose()*BB*U
print('singular quotient ranks',Aq.rank(),Bq.rank(),Aq.stack(Bq).rank(),(Aq+Bq).rank())

# Polynomial kernel over K(t), cleared to a primitive K[t]-vector.
Kt.<t> = PolynomialRing(K)
Ft=Kt.fraction_field()
Mf=matrix(Ft,Aq)+t*matrix(Ft,Bq)
print('generic quotient rank',Mf.rank())
ker=Mf.right_kernel_matrix()
print('generic quotient nullity',ker.nrows())
assert ker.nrows()==1
w=list(ker.row(0))
den=lcm([z.denominator() for z in w])
wp=[Kt(z*den) for z in w]
gcdw=gcd(wp)
wp=[z//gcdw for z in wp]
lead=next(z for z in wp if z)
wp=[z/lead.leading_coefficient() for z in wp]
print('kernel max degree',max(z.degree() for z in wp if z))
print('kernel support',sum(z!=0 for z in wp))
print('kernel quotient coordinates')
for j,z in enumerate(wp):
    if z: print(free[j],B6[comp[free[j]]],z)
assert (matrix(Kt,Aq)+t*matrix(Kt,Bq))*vector(Kt,wp)==0

# Lift to the original R_6 monomial basis.  This representative is enough
# because the common radical is killed by both A and B.
lift=[Kt(0)]*426
for j,z in enumerate(wp):
    lift[comp[free[j]]]=z
print('lifted monomial kernel')
for i,z in enumerate(lift):
    if z: print(B6[i],z)
assert (matrix(Kt,A)+t*matrix(Kt,B))*vector(Kt,lift)==0

# A cofactor in a coordinate where the primitive kernel is a nonzero constant
# detects every further rank drop of this symmetric corank-one pencil.
constant_coordinate=next(j for j,z in enumerate(wp) if z and z.degree()==0)
keep=[j for j in range(23) if j != constant_coordinate]
minor=(matrix(Kt,Aq)+t*matrix(Kt,Bq)).matrix_from_rows_and_columns(keep,keep)
det=minor.det()
print('constant kernel coordinate',constant_coordinate,wp[constant_coordinate])
print('22-principal-cofactor factor',det.factor())

# Do the same simultaneous-radical quotient on every sparsity component.
# Regular quotients are square, so their determinant gives the complete list
# of nonzero exceptional parameters.
print('all quotient determinants')
directed_ab_total=0
directed_ab_dim_total=0
directed_ba_total=0
directed_ba_dim_total=0
for comp0 in comps:
    if len(comp0)==1: continue
    A0=A.matrix_from_rows_and_columns(comp0,comp0)
    B0=B.matrix_from_rows_and_columns(comp0,comp0)
    N0=A0.stack(B0).right_kernel_matrix()
    piv0=N0.echelon_form().pivots()
    free0=[j for j in range(len(comp0)) if j not in piv0]
    U0=matrix(K,len(comp0),len(free0),lambda i,j: K(1) if i==free0[j] else K(0))
    A0q=U0.transpose()*A0*U0
    B0q=U0.transpose()*B0*U0
    KA0=A0q.right_kernel_matrix()
    KB0=B0q.right_kernel_matrix()
    directed_ab0=KA0*B0q*KA0.transpose()
    directed_ba0=KB0*A0q*KB0.transpose()
    directed_ab_total += directed_ab0.rank()
    directed_ab_dim_total += KA0.nrows()
    directed_ba_total += directed_ba0.rank()
    directed_ba_dim_total += KB0.nrows()
    if len(free0)==23:
        print('quotient',len(comp0),len(free0),'singular (handled above)',
              'directed A-to-B',directed_ab0.rank(),'of',KA0.nrows(),
              'directed B-to-A',directed_ba0.rank(),'of',KB0.nrows())
    else:
        det0=(matrix(Kt,A0q)+t*matrix(Kt,B0q)).det()
        print('quotient',len(comp0),len(free0),'ranks',A0q.rank(),B0q.rank(),
              'directed A-to-B',directed_ab0.rank(),'of',KA0.nrows(),
              'directed B-to-A',directed_ba0.rank(),'of',KB0.nrows(),
              'det',det0.factor())
print('directed A-to-B total',directed_ab_total,'of',directed_ab_dim_total)
print('directed B-to-A total',directed_ba_total,'of',directed_ba_dim_total)

# Factor the moving kernel as a degree-six polynomial.
X.<x0,x1,x2,x3,x4,x5> = PolynomialRing(K)
H0 = x3^2*(mu^2*x3^2 + 2*x3*x4 + (1-mu^2)*x4^2 + (mu^3-mu)*x4*x5)
H1 = x4*(x4^3 + mu^2*x3*x4^2 + (mu^3-mu)*x3^2*x5 + x3^3)
print('H0',H0)
print('H1',H1)
print('kernel factorization: (x2^2-a*x0*x1)*(H0+t*H1)')

v0=vector(Kt,[z[0] for z in lift])
v1=vector(Kt,[z[1] if z.degree()>=1 else 0 for z in lift])
assert matrix(Kt,A)*v0==0
assert matrix(Kt,B)*v1==0
assert matrix(Kt,A)*v1 + matrix(Kt,B)*v0==0
assert matrix(Kt,B)*v0 != 0
assert matrix(Kt,A)*v1 != 0
print('exact kernel chain verified: A*v0=0, A*v1+B*v0=0, B*v1=0')
