from itertools import product

# Clean exact calculation for the disjoint Aoki--Shioda pair on the Fermat
# sextic fourfold.  K=Q(mu,a), Phi_12(mu)=0, a^3=2.
C.<mu> = CyclotomicField(12)
Za.<aa> = PolynomialRing(C)
K.<a> = C.extension(aa^3-2)
mu=K(mu)
Kt.<t> = PolynomialRing(K)

TOP=(4,)*6
B6=[e for e in product(range(5),repeat=6) if sum(e)==6]
I6={e:i for i,e in enumerate(B6)}

def add_t(out,e,c):
    out[e]=out.get(e,Kt(0))+Kt(c)
    if out[e]==0: del out[e]

def mul_t(A,B):
    out={}
    for e,c in A.items():
        for f,d in B.items(): add_t(out,tuple(e[i]+f[i] for i in range(6)),c*d)
    return out

def add_k(out,e,c):
    out[e]=out.get(e,K(0))+K(c)
    if out[e]==0: del out[e]

def mul_k(A,B):
    out={}
    for e,c in A.items():
        for f,d in B.items(): add_k(out,tuple(e[i]+f[i] for i in range(6)),c*d)
    return out

def base_aoki():
    out={}; base=(1,0,0,2,0,0)
    fa=[((2,0,0,0,0,0),K(1)),((0,1,1,0,0,0),a)]
    fb=[((0,3,0,0,0,0),K(1)),((0,0,3,0,0,0),K(-1))]
    fc=[((0,0,0,0,j,4-j),mu^(4-j)) for j in range(5)]
    for (ea,ca),(eb,cb),(ec,cc) in product(fa,fb,fc):
        add_t(out,tuple(base[i]+ea[i]+eb[i]+ec[i] for i in range(6)),ca*cb*cc)
    return out

def permute_raw(poly,sigma):
    out={}
    for e,c in poly.items():
        ee=[0]*6
        for old,new in enumerate(sigma): ee[new]=e[old]
        add_t(out,tuple(ee),c)
    return out

def evaluate(poly,value):
    return {e:K(c(value)) for e,c in poly.items() if c(value)!=0}

def rho(poly):
    M=matrix(K,426,426,sparse=True)
    for i,left in enumerate(B6):
        for e,c in poly.items():
            right=tuple(TOP[j]-left[j]-e[j] for j in range(6))
            if min(right)>=0 and max(right)<=4 and right in I6:
                M[i,I6[right]]+=c
    return M

# Fermat Griffiths reduction from the indicated pole order to pole order 3.
def reduce_to_pole3(poly,pole):
    steps=pole-3
    denom=K(1)
    for current_pole in range(pole,3,-1):
        denom*=K(6*(current_pole-1))
    out={}
    for e,c in poly.items():
        counts=[e[i]//6 for i in range(6)]
        rem=[e[i]-6*counts[i] for i in range(6)]
        if sum(counts)!=steps or sum(rem)!=12 or any(v==5 for v in rem):
            continue
        numer=K(1)
        for i in range(6):
            for j in range(counts[i]): numer*=K(e[i]-5-6*j)
        add_k(out,tuple(rem),c*numer/denom)
    return out

def residue_pair(P,A):
    ans=K(0)
    for e,c in A.items():
        comp=tuple(TOP[i]-e[i] for i in range(6))
        if min(comp)>=0: ans+=c*P.get(comp,K(0))
    return ans

def equation_vector(P,poly,pole):
    vals=[]
    for q in B6:
        qpoly={tuple(q[i]+e[i] for i in range(6)):c for e,c in poly.items()}
        vals.append(residue_pair(P,reduce_to_pole3(qpoly,pole)))
    return vector(K,vals)

def compositions(total,parts):
    if parts==1:
        yield (total,); return
    for first in range(1,total-parts+2):
        for rest in compositions(total-first,parts-1):
            yield (first,)+rest

def solve_componentwise(M,rhs):
    adj=[set() for _ in B6]
    for (i,j),c in M.dict().items(): adj[i].add(j);adj[j].add(i)
    seen=set();sol=[K(0)]*426
    for seed in range(426):
        if seed in seen: continue
        seen.add(seed);todo=[seed];comp=[]
        while todo:
            u=todo.pop();comp.append(u)
            for v in adj[u]:
                if v not in seen: seen.add(v);todo.append(v)
        comp=sorted(comp)
        rc=vector(K,[rhs[i] for i in comp])
        if not any(rc): continue
        Mc=M.matrix_from_rows_and_columns(comp,comp)
        if Mc.augment(matrix(K,len(comp),1,list(rc))).rank()!=Mc.rank():
            return None
        hc=Mc.solve_right(rc)
        for j,idx in enumerate(comp): sol[idx]=hc[j]
    return {B6[i]:sol[i] for i in range(426) if sol[i]!=0}

# Raw P1 is obtained by substitution only.  sigma is odd, so the geometric
# residue of S' is -P1.  Thus raw t=+1 is [S]-[S'] and raw t=-1 is [S]+[S'].
P0=base_aoki()
P1=permute_raw(P0,(1,0,2,4,5,3))
Pt=dict(P0)
for e,c in P1.items(): add_t(Pt,e,t*c)

quad={(0,0,2,0,0,0):Kt(1),(1,1,0,0,0,0):-a}
H0={(0,0,0,4,0,0):mu^2,(0,0,0,3,1,0):2,
    (0,0,0,2,2,0):1-mu^2,(0,0,0,2,1,1):mu^3-mu}
H1={(0,0,0,0,4,0):1,(0,0,0,1,3,0):mu^2,
    (0,0,0,2,1,1):mu^3-mu,(0,0,0,3,1,0):1}
Ht={e:Kt(c) for e,c in H0.items()}
for e,c in H1.items(): add_t(Ht,e,t*c)
Gt=mul_t(quad,Ht)

# The self-pairing necessary obstruction <Pt,red(Gt^3)>.
G3t=mul_t(mul_t(Gt,Gt),Gt)
redt={}
for e,c in G3t.items():
    candidates=[i for i in range(6) if 6<=e[i]<=10
                and all(e[j]<=4 for j in range(6) if j!=i)]
    if len(candidates)==1:
        i=candidates[0];ee=list(e);ee[i]-=6
        add_t(redt,tuple(ee),c*K(e[i]-5)/K(18))
obstruction=Kt(0)
for e,c in redt.items():
    comp=tuple(TOP[i]-e[i] for i in range(6))
    if min(comp)>=0: obstruction+=c*Pt.get(comp,Kt(0))
expected=Kt(K(5)/K(2)*mu^2*a*t*(t+1))
assert obstruction==expected

print('RAW PENCIL (P0+t*P1_raw = [S]-t[S_prime])')
print('self_obstruction_factor =',obstruction.factor())

A=rho(evaluate(P0,0));B=rho(evaluate(P1,0))

def second_order_report(raw_t,label):
    P=evaluate(Pt,raw_t);G=evaluate(Gt,raw_t);M=A+K(raw_t)*B
    gvec=vector(K,[G.get(e,K(0)) for e in B6])
    assert M*gvec==0
    rhs0=equation_vector(P,mul_k(G,G),4)
    rhs=K(3)/K(2)*rhs0
    rank=M.rank();aug=M.augment(matrix(K,426,1,list(rhs))).rank()
    print(label)
    print('raw_t =',raw_t,'rho_rank =',rank,'second_augmented_rank =',aug)
    print('self_obstruction =',gvec*rhs0)
    return P,G,M,rhs

Pdiff,Gdiff,Mdiff,rhsdiff=second_order_report(1,'GEOMETRIC DIFFERENCE [S]-[S_prime]')
assert Mdiff.rank()==122
assert Mdiff.augment(matrix(K,426,1,list(rhsdiff))).rank()==123
assert vector(K,[Gdiff.get(e,K(0)) for e in B6])*equation_vector(Pdiff,mul_k(Gdiff,Gdiff),4)==5*mu^2*a

Psum,Gsum,Msum,rhssum=second_order_report(-1,'GEOMETRIC SUM [S]+[S_prime]')
assert Msum.rank()==122
assert Msum.augment(matrix(K,426,1,list(rhssum))).rank()==122

# Formal lifting of this chosen sum-class transverse vector through order 6.
hs={1:Gsum}
for order in range(2,7):
    total_rhs=vector(K,[0]*426)
    for power in range(2,order+1):
        coeff={}
        for js in compositions(order,power):
            term={(0,0,0,0,0,0):K(1)}
            for j in js: term=mul_k(term,hs[j])
            for e,c in term.items(): add_k(coeff,e,c)
        total_rhs+=K((-1)^power*(power+1))/K(2)*equation_vector(Psum,coeff,power+2)
    aug=Msum.augment(matrix(K,426,1,list(total_rhs))).rank()
    print('sum_formal_order =',order,'rho_rank =',Msum.rank(),'augmented_rank =',aug)
    assert aug==Msum.rank()
    hs[order]=solve_componentwise(Msum,total_rhs)
    assert hs[order] is not None
    hvec=vector(K,[hs[order].get(e,K(0)) for e in B6])
    assert Msum*hvec==total_rhs
    print('sum_correction_support =',len(hs[order]),'exact_residual_zero =',True)
