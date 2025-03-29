F = GF(17)
R.<x> = F[]

# P = R.random_element(degree=7)
P = 3*x^7 + 2*x^6 + 10*x^5 + 4*x^4 + 9*x^3 + 13*x^2 + 10*x + 2

# reed solomon extension by a factor of 4
DOMAIN_SIZE = P.degree() * 4

# split poly in even and odd coefficients
def split_polynomial(poly: R) -> tuple[R, R]:
    coeffs = poly.coefficients(sparse=False)
    deg = poly.degree()
    Pe = R(sum(coeffs[i] * x^(i//2) for i in range(0, deg + 1, 2)))
    Po = R(sum(coeffs[i] * x^((i-1)//2) for i in range(1, deg + 1, 2)))
    assert(Pe(x^2) + x * Po(x^2) == poly)
    return Pe, Po

# fold (commit) polynomial with random alpha
def fold_polynomial(poly: R, alpha: F) -> R:
    Pe, Po = split_polynomial(poly)
    Pf = Pe + alpha * Po
    return Pf

# evaluate poly over the domain
def evaluate_polynomial(poly: R) -> list[F]:
    return [poly(i) for i in range(DOMAIN_SIZE)]

# ROUND 1
evaluations1 = evaluate_polynomial(P)
# alpha should be random: `F.random_element()`
# but we'll use a fixed value for testing, for reproducibility
alpha1 = F(11)
folded1 = fold_polynomial(P, alpha1)

# ROUND 2
evaluations2 = evaluate_polynomial(folded1)
alpha2 = F(6)
folded2 = fold_polynomial(folded1, alpha2)

# ROUND 3
evaluations3 = evaluate_polynomial(folded2)
alpha3 = F(13)
folded3 = fold_polynomial(folded2, alpha3)

# we stop here because the degree is 0 (a constant)
assert(folded3.degree() == 0)

# Now to the verification

# First we pick a random index z
# again I hardcoded it for the example
z = F(8)

# get P(z) and P(-z)
P_z = evaluations1[z]
P_z_minus = evaluations1[-z]

# in the next layers, we square the index at each step
z2 = z^2
P_z2 = evaluations2[z2]
P_z_minus2 = evaluations2[-z2]

z4 = z2^2
P_z4 = evaluations3[z4]
P_z_minus4 = evaluations3[-z4]

round1_check = (((alpha1 + z) / (2*z)) * P_z) + (((alpha1 - z) / (2*(-z))) * P_z_minus)
assert(round1_check == P_z2)

round2_check = (((alpha2 + z2) / (2*z2)) * P_z2) + (((alpha2 - z2) / (2*(-z2))) * P_z_minus2)
assert(round2_check == P_z4)

round3_check = (((alpha3 + z4) / (2*z4)) * P_z4) + (((alpha3 - z4) / (2*(-z4))) * P_z_minus4)

# in the last check, we verify that it's equal to the last folded polynomial
# a constant, which is known by the verifier
assert(round3_check == folded3)

print("VERIFIED!")