
GENERATOR = 2

def gf_mult(x, y, prim=0, field_charac_full=256):
    res = 0

    while y:
        if y & 1==1:
            res = res ^ x
        y = y >> 1  # y // 2
        x = x << 1  # x * 2

        if prim > 0 and x & field_charac_full:  # if x ~ 256 use xor
            x = x ^ prim

    return res


def find_prime_polys(c_exp=8):                                             # ищет неприводимые полиномы
    root_charac = 2
    field_charac = pow(root_charac, c_exp) - 1  # 255
    field_charac_next = pow(root_charac, c_exp + 1) - 1  # 511

    prim_candidates = []
    irreducible_polynomial = []
    seen = [0] * (field_charac + 1)

    for i in range(field_charac + 2, field_charac_next, root_charac):  # we need the number > 256 to avoid repetitions
        prim_candidates.append(i)  # try every simple polynomial, excluding even ones

    # Here is implemented a bruteforce approach to find all these prime polynomials, by generating every possible
    # prime polynomials(ie, every integers between field_charac + 1 and field_charac * 2), and then we build
    # the whole Galois Field, and we reject the candidate prime polynomial if it duplicates even one value or if it
    # generates a value above field_charac(ie, cause an overflow).
    for i in range(len(prim_candidates)):
        prim = prim_candidates[i]

        seen = [0] * (field_charac + 1)  # specifies whether the value in the field has already been generated (seen[x] ?== 1) or not

        conflict = False  # to find out if there was at least 1 conflict
        x = 1

        for j in range(field_charac):
            x = gf_mult(x, GENERATOR, prim, field_charac + 1)  # compute the next value in the field

            if x > field_charac or seen[x] == 1:  # if this number is a duplicate, then we reject it
                conflict = True
                break
            else:
                seen[x] = 1  # remember this value to detect future duplicates

        if not conflict:  # if there was no conflict, then this is a simple polynomial
            irreducible_polynomial.append(prim)

    return irreducible_polynomial




def init_tables(prim=285, c_exp=8):                     # получим таблицы
    field_charac = 2 ** c_exp - 1  # 255

    gf_exp = [0] * (field_charac * 2)
    gf_log = [0] * (field_charac + 1)

    # For each element from the Galois field, we calculate the log and exp
    x = 1
    for i in range(field_charac):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult(x, GENERATOR, prim, field_charac + 1)

    for i in range(field_charac, field_charac * 2):
        gf_exp[i] = gf_exp[i - field_charac]

    return gf_exp, gf_log















