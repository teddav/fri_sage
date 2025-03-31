F = GF(0x10001)
R.<x> = F[]

# x == 0 for x in [0,...,9]
Z = prod((x - i) for i in range(10))

# first we make sure P == 0 for x < 10
points = [[i, 0] for i in range(10)]
# then we just add random points up to 100
points += [[i, F.random_element()] for i in range(10, 101)]
P = R.lagrange_polynomial(points)

assert(P.degree() == 100)
assert(P % Z == 0)

Q = P // Z
assert(Q.degree() == P.degree() - Z.degree())

# the evaluation domain is normally based on the root of unity
EVALUATION_DOMAIN = range(1000, 2000)
def evaluate(poly):
    return [poly(i) for i in EVALUATION_DOMAIN]

evaluations_P = evaluate(P)
evaluations_Z = evaluate(Z)
evaluations_Q = evaluate(Q)

# we check at 5 different points
for _ in range(5):
    z = Integer(F.random_element()) % 1000
    assert(evaluations_P[z] == evaluations_Q[z] * evaluations_Z[z])

# we generate a random P
P_fake = R.random_element(degree=100)
evaluations_P_fake = evaluate(P_fake)
evaluations_Q_fake = [P_fake(i) / Z(i) for i in EVALUATION_DOMAIN]

for _ in range(5):
    z = Integer(F.random_element()) % 1000
    assert(evaluations_P_fake[z] == evaluations_Q_fake[z] * evaluations_Z[z])

# we recover Q from the fake evaluations
Q_fake = R.lagrange_polynomial(zip(EVALUATION_DOMAIN, evaluations_Q_fake))
print("Q' degree = ", Q_fake.degree())
