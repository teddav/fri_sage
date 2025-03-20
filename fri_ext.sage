# docker run --rm --platform linux/amd64 -v $(pwd):/app -w /app sagemath/sagemath 'sage ./fri_article.sage'

from merkle import MerkleTree, hash

F = GF(97)
R.<u> = F[]

# quadratic extension over F
M = R(u^2 + 96*u + 5) # irreducible polynomial: R.irreducible_element(2)
EXT.<v> = F.extension(M)

# we need a way to convert numbers in our extension to byte
# this is needed in order to hash them for the Merkle tree
def ext_to_bytes(v: EXT) -> bytes:
    return "|".join([str(c) for c in v.polynomial().coefficients()]).encode()

OMEGA = F.multiplicative_generator()

P = R.random_element(degree=9)
print("P", P)

P_degree = P.degree()
BLOWUP_FACTOR = 4

# split poly in even and odd coefficients
def split_polynomial(poly: R) -> tuple[R, R]:
    x = poly.parent().gen()
    coeffs = poly.coefficients(sparse=False)
    deg = poly.degree()
    Pe = sum(coeffs[i] * x^(i//2) for i in range(0, deg + 1, 2))
    Po = sum(coeffs[i] * x^((i-1)//2) for i in range(1, deg + 1, 2))
    assert(Pe(x^2) + x * Po(x^2) == poly)
    return Pe, Po

# fold (commit) polynomial with random alpha
def fold_polynomial(poly: R, alpha: EXT) -> R:
    Pe, Po = split_polynomial(poly)
    Pf = Pe + alpha * Po
    return Pf

# evaluate poly over the domain
def evaluate_polynomial(poly: R, blowup: int, omega: F) -> list[F]:
    domain = (P_degree + 1) * blowup
    return [poly(omega^i) for i in range(domain)]

# commitment: compute the merkle tree of the evaluations
def commit_polynomial(poly: R, blowup: int, omega: F) -> MerkleTree:
    evaluations = evaluate_polynomial(poly, blowup, omega)
    leaves = [hash(ext_to_bytes(eval)) for eval in evaluations]
    return MerkleTree(leaves)

# one round of FRI commitment:
# - compute the merkle tree
# - fold the polynomial
def round(poly: R, omega: F) -> tuple[MerkleTree, R, EXT]:
    tree = commit_polynomial(poly, BLOWUP_FACTOR, omega)
    alpha = EXT.random_element()
    folded = fold_polynomial(poly, alpha)
    return (tree, folded, alpha)

# fold and commit the poly until we reach a constant (degree 0):
def commit(poly: R) -> tuple[list[MerkleTree], list[R], list[EXT]]:
    trees = []
    polys = []
    alphas = []
    omega = OMEGA
    
    polys.append(poly)

    while poly.degree() > 0:
        (tree, folded, alpha) = round(poly, omega)
        trees.append(tree)
        alphas.append(alpha)
        polys.append(folded)

        poly = folded
        omega = pow(omega, 2)

    return (trees, polys, alphas)

# open commitment at index `z`
# that means:
# - computing f(z) and f(-z)
# - generating the merkle proof
def query(z: F, polys: list[R], trees: list[MerkleTree]) -> list[tuple[EXT, EXT, list[bytes]]]:
    omega = OMEGA ^ z
    queries = []
    for i in range(len(trees)):
        f_z = polys[i](omega)
        f_z_minus = polys[i](-omega)

        merkle_proof = trees[i].proof(z)
        queries.append((f_z, f_z_minus, merkle_proof))
        omega = pow(omega, 2)
    return queries

# verify that all the layers match
# refer to the article if you don't understand the formula
def verify(z: F, queries: list[tuple[EXT, EXT, list[bytes]]], commitments: list[bytes], alphas: list[EXT]):
    omega = OMEGA ^ z
    previous_round = None
    for i in range(len(queries)):
        f_z = queries[i][0]

        merkle_proof = queries[i][2]
        root = MerkleTree.compute_tree_from_proof(ext_to_bytes(f_z), merkle_proof)
        assert root == commitments[i]

        f_z_minus = queries[i][1]
        verif = (((alphas[i] + omega) / (2*omega)) * f_z) + (((alphas[i] - omega) / (2*(-omega))) * f_z_minus)
        
        if previous_round:
            assert f_z == previous_round
        
        previous_round = verif
        omega = pow(omega, 2)
    
    # last round
    assert previous_round == commitments[-1]


(trees, polys, alphas) = commit(P)
commitments = [tree.root() for tree in trees] + [polys[-1]]

z = 12 # random
assert(z <= BLOWUP_FACTOR * P_degree)

queries = query(z, polys, trees)
verify(z, queries, commitments, alphas)