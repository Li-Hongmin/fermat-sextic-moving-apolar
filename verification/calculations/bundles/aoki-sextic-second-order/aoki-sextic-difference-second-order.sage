from itertools import product

# Exact low-memory certificate for the mixed second-order obstruction of the
# disjoint Aoki--Shioda pair on the Fermat sextic fourfold.  All linear algebra
# is split into the simultaneous sparse blocks of A, B, C_G, and C_kappa.
C0.<mu0> = CyclotomicField(12)
Za.<aa> = PolynomialRing(C0)
K.<a> = C0.extension(aa^3 - 2)
mu = K(mu0)

TOP = (4,)*6
B6 = [e for e in product(range(5), repeat=6) if sum(e) == 6]
I6 = {e:i for i,e in enumerate(B6)}

def addterm(out, e, c):
    c = K(c)
    out[e] = out.get(e, K(0)) + c
    if out[e] == 0:
        del out[e]

def mul(left, right):
    out = {}
    for e,c in left.items():
        for f,d in right.items():
            addterm(out, tuple(e[i]+f[i] for i in range(6)), c*d)
    return out

def base_aoki():
    out = {}
    base = (1,0,0,2,0,0)
    fa = [((2,0,0,0,0,0),K(1)), ((0,1,1,0,0,0),a)]
    fb = [((0,3,0,0,0,0),K(1)), ((0,0,3,0,0,0),K(-1))]
    fc = [((0,0,0,0,j,4-j),mu^(4-j)) for j in range(5)]
    for (ea,ca),(eb,cb),(ec,cc) in product(fa,fb,fc):
        exponent = tuple(base[i]+ea[i]+eb[i]+ec[i] for i in range(6))
        addterm(out, exponent, ca*cb*cc)
    return out

def permute_variables(poly, sigma):
    out = {}
    for e,c in poly.items():
        exponent = [0]*6
        for old,new in enumerate(sigma):
            exponent[new] = e[old]
        addterm(out, tuple(exponent), c)
    return out

def rho(poly):
    out = matrix(K, len(B6), len(B6), sparse=True)
    for row,left in enumerate(B6):
        for e,c in poly.items():
            right = tuple(TOP[j]-left[j]-e[j] for j in range(6))
            col = I6.get(right)
            if col is not None:
                out[row,col] += c
    return out

# C_H(q,w)=<P,R_{4->3}(q H w)> for H in degree six.  For a monomial p
# of P, reduction in x_i forces w=TOP-p-q-h+6e_i and contributes
# (5-p_i)/18.  This directly emits the sparse matrix, without polynomial
# expansion in degree 24.
def mixed_matrix(P, H):
    out = matrix(K, len(B6), len(B6), sparse=True)
    inv18 = K(18)^(-1)
    for p,cp in P.items():
        for h,ch in H.items():
            target0 = [TOP[j]-p[j]-h[j] for j in range(6)]
            for var in range(6):
                coefficient = cp*ch*K(5-p[var])*inv18
                if coefficient == 0:
                    continue
                target = list(target0)
                target[var] += 6
                for row,q in enumerate(B6):
                    w = tuple(target[j]-q[j] for j in range(6))
                    col = I6.get(w)
                    if col is not None:
                        out[row,col] += coefficient
    return out

def simultaneous_blocks(matrices):
    adjacency = [set() for _ in B6]
    for mat in matrices:
        for (i,j),coefficient in mat.dict().items():
            if coefficient != 0:
                adjacency[i].add(j)
                adjacency[j].add(i)
    seen = set()
    blocks = []
    for seed in range(len(B6)):
        if seed in seen:
            continue
        seen.add(seed)
        todo = [seed]
        block = []
        while todo:
            i = todo.pop()
            block.append(i)
            for j in adjacency[i]:
                if j not in seen:
                    seen.add(j)
                    todo.append(j)
        blocks.append(sorted(block))
    return sorted(blocks, key=lambda block:(-len(block),block))

P0 = base_aoki()
P1 = permute_variables(P0, (1,0,2,4,5,3))
P = dict(P0)
for e,c in P1.items():
    addterm(P, e, c)

quad = {(0,0,2,0,0,0):K(1), (1,1,0,0,0,0):-a}
H = {
    (0,0,0,4,0,0):mu^2,
    (0,0,0,3,1,0):K(3),
    (0,0,0,2,2,0):K(1)-mu^2,
    (0,0,0,2,1,1):K(2)*(mu^3-mu),
    (0,0,0,0,4,0):K(1),
    (0,0,0,1,3,0):mu^2,
}
G = mul(quad, H)

# The compact strong certificate found in the unique excess block:
# kappa=x2^2(x5^4+mu^2*x4^4+(1-mu^2)*x3^4).
kappa = {
    (0,0,2,0,0,4):K(1),
    (0,0,2,0,4,0):mu^2,
    (0,0,2,4,0,0):K(1)-mu^2,
}

A = rho(P0)
B = rho(P1)
CG = mixed_matrix(P, G)
Ck = mixed_matrix(P, kappa)
assert CG == CG.transpose()
assert Ck == Ck.transpose()

gvec = vector(K, [G.get(e,K(0)) for e in B6])
kvec = vector(K, [kappa.get(e,K(0)) for e in B6])

blocks = simultaneous_blocks((A,B,CG,Ck))
max_block = max(map(len,blocks))
common_kernel_dimension = 0
mixed_support = 0
pure_support = 0
rank_A = 0
rank_B = 0
rank_stack = 0
rank_sum = 0

# A, B, CG, and Ck are simultaneously block diagonal for these blocks.
# Construct common-kernel bases only inside one block at a time.
for block in blocks:
    Ab = A.matrix_from_rows_and_columns(block,block)
    Bb = B.matrix_from_rows_and_columns(block,block)
    CGb = CG.matrix_from_rows_and_columns(block,block)
    Ckb = Ck.matrix_from_rows_and_columns(block,block)
    Wb = Ab.stack(Bb).right_kernel_matrix()
    Mb = Ab+Bb
    common_kernel_dimension += Wb.nrows()
    rank_A += Ab.rank()
    rank_B += Bb.rank()
    rank_stack += len(block)-Wb.nrows()
    rank_sum += Mb.rank()

    gb = vector(K, [gvec[i] for i in block])
    kb = vector(K, [kvec[i] for i in block])
    assert Mb*gb == 0
    assert Mb*kb == 0
    mixed = Wb*(kb*CGb).column()
    mixed_support += sum(z != 0 for z in mixed)
    assert mixed == 0

    # Stream the quadratic restriction one common-kernel vector at a time;
    # never materialize the global 303-by-303 restriction (nor even a dense
    # block-sized restriction).
    for w in Wb.rows():
        pure_column = Wb*(Ckb*w.column())
        pure_support += sum(z != 0 for z in pure_column)
        assert pure_column == 0

pairing = kvec*CG*gvec
expected = K(2)/K(3)*(K(1)-mu^2)*a
assert pairing == expected

print('field_relations: mu^4-mu^2+1=0, a^3=2')
print('kappa: x2^2*(x5^4+mu^2*x4^4+(1-mu^2)*x3^4)')
print('supports(P,G,kappa):', len(P), len(G), len(kappa))
print('ranks(A,B,[A;B],A+B):', rank_A, rank_B, rank_stack, rank_sum)
print('joint_blocks/count/max:', len(blocks), max_block)
print('common_kernel_dimension:', common_kernel_dimension)
print('M*kappa:', True)
print('kappa*C_G*ker(A;B):', mixed_support)
print('ker(A;B)*C_kappa*ker(A;B)^T:', pure_support)
print('kappa*C_G*G:', pairing)
