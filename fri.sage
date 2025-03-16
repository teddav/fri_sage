# docker run --rm --platform linux/amd64 -v $(pwd):/app -w /app sagemath/sagemath 'sage ./fri.sage'

from merkle import MerkleTree, hash

F = GF(97)
R.<x> = F[]
OMEGA = F.multiplicative_generator()
P = R.random_element(degree=9)
print("P", P)
P_degree = P.degree()

BLOWUP_FACTOR = 4

def split_polynomial(poly: R) -> tuple[R, R]:
    x = poly.parent().gen()
    coeffs = poly.coefficients(sparse=False)
    deg = poly.degree()
    Pe = R(sum(coeffs[i] * x^(i//2) for i in range(0, deg + 1, 2)))
    Po = R(sum(coeffs[i] * x^((i-1)//2) for i in range(1, deg + 1, 2)))
    return Pe, Po

def fold_polynomial(poly: R, alpha: F) -> R:
    Pe, Po = split_polynomial(poly)
    Pf = Pe + alpha * Po
    return Pf

def evaluate_polynomial(poly: R, blowup: int, omega: F) -> list[F]:
    # size = (poly.degree() + 1) * blowup
    size = (P_degree + 1) * blowup
    return [poly(omega^i) for i in range(size)]

def commit_polynomial(poly: R, blowup: int, omega: F) -> MerkleTree:
    evaluations = evaluate_polynomial(poly, blowup, omega)
    leaves = [hash(eval.to_bytes()) for eval in evaluations]
    return MerkleTree(leaves)

def round(poly: R, omega: F) -> tuple[MerkleTree, R, F]:
    tree = commit_polynomial(poly, BLOWUP_FACTOR, omega)
    r = F.random_element()
    folded = fold_polynomial(poly, r)
    return (tree, folded, r)

def commit(poly: R) -> tuple[list[MerkleTree], list[R], list[F]]:
    trees = []
    polys = []
    randoms = []
    omega = OMEGA
    while poly.degree() > 0:
        polys.append(poly)
        (tree, folded, r) = round(poly, omega)
        trees.append(tree)
        randoms.append(r)
        poly = folded
        omega = pow(omega, 2)
    polys.append(poly)

    # last round
    # (tree, folded, r) = round(poly, omega)
    # trees.append(tree)
    # polys.append(folded)
    # randoms.append(r)

    return (trees, polys, randoms)

def query(index: F, polys: list[R], trees: list[MerkleTree]) -> list[tuple[F, F, list[bytes]]]:
    omega = OMEGA ^ index
    queries = []
    for i in range(len(trees)):
        # if index >= (polys[i].degree() + 1) :
        #     index = index % (polys[i].degree() + 1)
        f_g = polys[i](omega)
        f_g_minus = polys[i](-omega)
        merkle_proof = trees[i].proof(index)
        queries.append((f_g, f_g_minus, merkle_proof))
        omega = pow(omega, 2)
    return queries

def verify(index: F, queries: list[tuple[F, F, list[bytes]]], commitments: list[bytes], randoms: list[F]):
    omega = OMEGA ^ index
    previous_round = None
    for i in range(len(queries)):
        f_g = queries[i][0]

        merkle_proof = queries[i][2]
        root = MerkleTree.compute_tree_from_proof(f_g.to_bytes(), merkle_proof)
        assert root == commitments[i]

        f_g_minus = queries[i][1]
        verif = (((randoms[i] + omega) / (2*omega)) * f_g) + (((randoms[i] - omega) / (2*(-omega))) * f_g_minus)
        
        if previous_round:
            assert f_g == previous_round
        
        previous_round = verif
        omega = pow(omega, 2)
    
    # last round
    assert previous_round == commitments[-1]


(trees, polys, randoms) = commit(P)
print("COMMIT")
commitments = [tree.root() for tree in trees] + [polys[-1]]
index = 30
assert(index <= BLOWUP_FACTOR * P_degree)
print("QUERY")
queries = query(index, polys, trees)
print("VERIFY")
verify(index, queries, commitments, randoms)
