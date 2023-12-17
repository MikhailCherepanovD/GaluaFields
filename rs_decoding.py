from gf import *
from help_gf import *


''' @brief Function calculates the polynomial of error syndromes
	* @param msg - an incoming message that is represented by a vector of polynomials(in this case, integers)
	* @param red_code_len - the number of characters representing the redundant code
	* @return generator polynomial '''

def rs_calc_syndromes(msg, red_code_len):                     # вычисляет синдромы ошибок
    synd = [0] * (red_code_len + 1)

    for i in range(1, red_code_len + 1):
        temp = pow_my(GENERATOR, i - 1)
        synd[i] = gf_poly_eval(msg, temp)

    return synd




''' @brief Function find the error locator polynomial L(x) according to the known location of the error (for correct errarta)
	* @param err_pos - vector with the positions of the erroneous characters
	* @return error locator polynomial L(x) 
	'''


def rs_find_errarta_locator(err_pos):                     # находит локаторы ошибок
    e_loc = [1]
    temp = [1]
    temp2 = [0, 0]
    add = []

    for i in range(len(err_pos)):
        temp2[0] = pow_my(GENERATOR, err_pos[i])  # since we know location of the error, we can find L(x) as
        add = gf_poly_add(temp, temp2)  #
        e_loc = gf_poly_mult(e_loc, add)

    return e_loc



''' @brief Function find the error polynomial W(x) = S(x)*L(x)
	* @param synd - polynomial of error syndromes (vector int) S(x)
	* @param err_loc - error locator polynomial L(x)
	* @param err_loc_size - size of L(x)
	* @return error polynomial W(x) 
	'''
def rs_find_error_evaluator(synd, err_loc, err_loc_size):
    poly_mul = gf_poly_mult(synd, err_loc)
    remainder = poly_mul[-err_loc_size:]  # переменная remainder будет содержать последние err_loc_size элементов из списка poly_mul.

    # Since W(x) cannot exceed u-1, where u is the number of errors, we will use division to discard the extra part
    # gf_poly_div(poly_mul, help, quotient, remainder)

    return remainder








def rs_correct_errata(msg_in, synd, err_pos):
    coef_pos = [len(msg_in) - 1 - pos for pos in err_pos]

    # find the error locator polynomial L(x) according to the known location of the error
    err_loc = rs_find_errarta_locator(coef_pos)

    # find the error polynomial W(x)
    synd.reverse()
    err_eval = rs_find_error_evaluator(synd, err_loc, len(err_loc))
    err_eval.reverse()

    # x - will store the position of the errors
    # we need to get the error location polynomial X from the error positions in err_pos
    # (the roots of the error locator polynomial, i.e., where it evaluates to 0)
    x = [pow_my(GENERATOR, -(255 - coef_pos[i])) for i in range(len(coef_pos))]

    E = [0] * len(msg_in)  # will store the values that need to be corrected to the original message with errors

    Xlength = len(x)
    err_eval.reverse()



    for i in range(len(x)):
        Xi_inv = inverse(x[i])

        # Find the formal derivative of the error locator polynomial
        # the formal derivative of the errata locator is used as the denominator of the Forney Algorithm,
        # which simply says that the ith error value is given by error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)).
        err_loc_prime_tmp = [sub(1, mult(Xi_inv, x[j])) for j in range(Xlength) if j != i]

        # error polynomial Yi = W(Xi^(-1))/L'(Xi^(-1))
        # calculate the product that is the denominator of the Forney algorithm (the derivative of the error locator)
        err_loc_prime = 1
        for coef in err_loc_prime_tmp:
            err_loc_prime = mult(err_loc_prime, coef)

        err_loc_prime_tmp.clear()

        y = gf_poly_eval(err_eval, Xi_inv)  # numerator
        y = mult(pow_my(x[i], 1), y)

        try:
            if err_loc_prime == 0:  # divisor should not be 0
                raise ValueError("Could not find error magnitude")
        except ValueError as ex:
            print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
            print("\nError:", ex)
            print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n")
            exit(-1)

        magnitude = div_my(y, err_loc_prime)  # The value of the error value calculated by the Forney algorithm
        # Dividing the error estimator by the derivative of the error locator
        E[err_pos[i]] = magnitude  # gives us the error value, that is, the value for recovering the symbol



    msg_in = gf_poly_add(msg_in, E);                                   # C(x) = C'(x) + E(x) (xor)
    return msg_in






def rs_find_error_locator(synd, red_code_len):
    err_loc = [1]  # error locator polynomial (sigma) C(x)
    old_loc = [1]  # the error locator polynomial of the previous iteration

    synd_shift = len(synd) - red_code_len

    # The Berlekamp–Massey algorithm is an alternative to the Reed–Solomon Peterson decoder for solving the set of linear equations.
    # The main idea is that the algorithm iteratively evaluates the error locator polynomial. To do this, it calculates the delta discrepancy,
    # by which we can determine whether we need to update the error locator or not
    k, delta = 0, 0

    for i in range(red_code_len):
        k = i + synd_shift

        # calculating the delta of the discrepancy
        delta = synd[k]

        for j in range(1, len(err_loc)):                                                    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            delta ^= mult(err_loc[len(err_loc) - 1 - j], synd[k - j])  # delta = Sn + C1*Sn-1 +..+ Cj*Sk-j

        # shift the polynomials to calculate the next degree
        old_loc.append(0)


        if delta != 0:  # if delta == 0, the algorithm assumes that C(x) and L are correct for the moment and continues
            if len(old_loc) > len(err_loc):  # ~2*L <= k + erase_count
                # Computing errata locator polynomial Sigma
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, inverse(delta))
                err_loc = new_loc

            # Update with the discrepancy
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))                                    #!!!!!!!!!!!!!!!!!!!!



    while err_loc and err_loc[0] == 0:
        err_loc.pop(0)

    errs = len(err_loc) - 1
    try:
        if errs * 2 > red_code_len:
            raise ValueError("Too many errors to correct")
    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("\nError:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        exit(-1)

    return err_loc


def rs_find_errors(err_loc, nmess):
    err_pos = []

    errs = len(err_loc) - 1
    for i in range(nmess):
        if gf_poly_eval(err_loc, pow_my(GENERATOR, i)) == 0:
            err_pos.append(nmess - 1 - i)

    try:
        if len(err_pos) != errs:
            raise ValueError("Too many (or few) errors found for the errata locator polynomial!")
    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("\nError:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n")
        exit(-1)

    return err_pos



def rs_decode_msg(msg_in, red_code_len):
    try:
        if len(msg_in) > 255:
            raise ValueError("Message is too long")
    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("\nError:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n")
        exit(-1)

    msg_out = msg_in

    # so that we do not count the generator polynomial several times and do not divide,
    # we immediately count the error syndrome polynomial, and if there is not at least
    # one non-0 value in it, then the message is not distorted
    synd = rs_calc_syndromes(msg_out, red_code_len)
    max_val = max(synd)

    if max_val == 0:
        return msg_out




    # Find the error locator polynomial L(x)
    err_loc = rs_find_error_locator(synd, red_code_len)

    err_loc.reverse()

    # Find the vector of the index of the characters that need to be corrected
    err_pos = rs_find_errors(err_loc, len(msg_out))

    try:
        if not err_pos:
            raise ValueError("Could not locate error")
    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("\nError:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n")
        exit(-1)



    # Find errors values and apply them to correct the message
    # Compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
    msg_out = rs_correct_errata(msg_out, synd, err_pos)

    # Count the error syndrome polynomial, and if there is not at least
    # one non-0 value in it, then the message is decoded successfully
    synd = rs_calc_syndromes(msg_out, red_code_len)
    max_synd = max(synd)
    try:
        if max_synd > 0:
            raise ValueError("Could not correct message")
    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("\nError:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n")
        exit(-1)

    return msg_out


















