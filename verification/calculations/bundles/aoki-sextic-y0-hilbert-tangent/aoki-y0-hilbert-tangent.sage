# First-order Hilbert tangent probe for the central Y0 ideal.
from itertools import product
import os
if os.environ.get('Y0_EXACT')=='1':
    C.<mu0>=CyclotomicField(12)
    Ta.<aa>=PolynomialRing(C)
    K.<a0>=C.extension(aa^3-2)
    F=K.absolute_field('t')
    toF=F.structure()[1]
    mu=toF(mu0); a=toF(a0); ic=mu^3
else:
    p=109; F=GF(p); mu=F(8); a=F(57); ic=F(76)
R.<x0,x1,x2,x3,x4,x5>=PolynomialRing(F)
L0=x4-mu*x5; L1=x5-mu*x3
Q0=x0^2-a*x1*x2; Q1=x1^2-a*x0*x2
C0=x1^3+x2^3+ic*x3^3; C1=x0^3+x2^3+ic*x4^3
F0=sum(x^6 for x in R.gens())
I=R.ideal([L0*L1,Q0*L1,Q1*L0,C1*L0,C0*L1,F0])
G=I.groebner_basis(); print('Y0 ideal groebner length',len(G))

def basis(d):
    return [e for e in product(range(d+1),repeat=6) if sum(e)==d]
def mon(e): return prod(R.gen(i)^e[i] for i in range(6))

degs=[2,3,3,4,4,6]
Syz=I.syzygy_module(); print('syzygy shape',Syz.nrows(),Syz.ncols())
syz_rows=[list(Syz.row(i)) for i in range(Syz.nrows())]
unknown=[]; posmap={}
for i,d in enumerate(degs):
    for j,e in enumerate(basis(d)):
        unknown.append((i,j,e)); posmap[(i,j)]=len(unknown)-1
U=len(unknown); eqdict={}
for rel in syz_rows:
    rdeg=None
    for i,c in enumerate(rel):
        if c: rdeg=degs[i]+c.degree(); break
    if rdeg is None: continue
    for i,c in enumerate(rel):
        if not c: continue
        for j,e in enumerate(basis(degs[i])):
            term=(c*mon(e)).reduce(G); col=posmap[(i,j)]
            for target_e,val in term.dict().items():
                row=eqdict.setdefault((rdeg,target_e),{})
                row[col]=row.get(col,F(0))+val
rows=[]
for sparse_row in eqdict.values():
    row=[F(0)]*U
    for col,val in sparse_row.items(): row[col]=val
    if any(row): rows.append(row)
M=matrix(F,rows,sparse=True)
matrix_rank=M.rank()
print('Polynomial representative ambient dim',U,'equations',M.nrows(),'rank',matrix_rank,'kernel dim',U-matrix_rank)
leading=[tuple(g.lm().exponents()[0]) for g in G]
def ideal_piece_dim(d):
    return sum(any(all(e[j]>=l[j] for j in range(6)) for l in leading) for e in basis(d))
zero_dims=[ideal_piece_dim(d) for d in degs]
assert zero_dims==[1,8,8,33,33,225]
print('Zero representative dimensions',zero_dims,'total',sum(zero_dims))
print('Embedded Hilbert tangent dimension',U-matrix_rank-sum(zero_dims))

h1=(-x2^2*x4^4-mu^2*x2^2*x3*x4^3+(1-mu^2)*x2^2*x3^2*x4^2
    +x2^2*x3^3*x4+mu^2*x2^2*x3^4+a*x0*x1*x4^4+mu^2*a*x0*x1*x3*x4^3
    +(mu^2-1)*a*x0*x1*x3^2*x4^2-a*x0*x1*x3^3*x4-mu^2*a*x0*x1*x3^4)
extra_rows=[]; target=[]; hr=h1.reduce(G)
for target_e in basis(6):
    row=[F(0)]*U
    for j,e in enumerate(basis(6)):
        term=mon(e).reduce(G)
        row[posmap[(5,j)]]=term.monomial_coefficient(mon(target_e))
    tv=hr.monomial_coefficient(mon(target_e))
    if any(row): extra_rows.append(row); target.append(tv)
    else: assert tv==0
Minc=matrix(F,rows+extra_rows,sparse=True)
binc=vector(F,[0]*len(rows)+target)
print('incidence h1 rank/aug',Minc.rank(),Minc.augment(matrix(F,len(binc),1,list(binc))).rank(),flush=True)
