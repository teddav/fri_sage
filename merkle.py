from hashlib import sha256

def hash(data: bytes) -> bytes:
        return sha256(data).digest()

def hash_leaves(left: bytes, right: bytes) -> bytes:
        return hash(left + right)

class MerkleTree:
    def __init__(self, leaves: list[bytes]):
        self.tree = self.compute_tree(leaves)

    def root(self) -> bytes:
        return self.tree[-1][0]

    def compute_tree(self, leaves: list[bytes]) -> bytes:
        if len(leaves) % 2 != 0:
            leaves.append(hash(b"0"))
        tree = [leaves]

        while len(leaves) != 1:
            leaves = [hash_leaves(leaves[i], leaves[i+1]) for i in range(0, len(leaves), 2)]
            if len(leaves) % 2 != 0 and len(leaves) != 1:
                leaves.append(hash(b"0"))
            tree.append(leaves)
        return tree
    
    def proof(self, leaf_index: int) -> bytes:
        proof = []
        for i in range(len(self.tree) - 1):
            if leaf_index % 2 == 0:
                proof.append((0, self.tree[i][leaf_index + 1]))
            else:
                proof.append((1, self.tree[i][leaf_index - 1]))
            leaf_index //= 2
        return proof
    
    def compute_tree_from_proof(leaf: bytes, proof: list[bytes]) -> bytes:
        leaf = hash(leaf)

        for (direction, sibling) in proof:
            if direction == 0:
                leaf = hash_leaves(leaf, sibling)
            else:
                leaf = hash_leaves(sibling, leaf)
        return leaf

def test():
    leaves = [
        b"123",
        b"456",
        b"789",
        b"1011",
    ]
    hashed = [hash(leaf) for leaf in leaves]
    mt = MerkleTree(hashed)
    print(mt.root().hex())
    proof = mt.proof(0)
    print(proof)
    print(mt.compute_tree_from_proof(leaves[0], proof).hex())
