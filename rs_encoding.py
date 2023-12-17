from gf import *
from help_gf import *



def rs_generator_poly(red_code_len):
    g = [1]
    temp = [1, 0]

    for i in range(red_code_len):
        temp[1] = pow_my(GENERATOR, i)
        g = gf_poly_mult(g, temp)

    return g




def rs_encode_msg(msg_in, red_code_len):
    try:
        if len(msg_in) + red_code_len < 256:
            msg_in_size = len(msg_in)
            gen = rs_generator_poly(red_code_len)
            msg_out_size = msg_in_size + red_code_len

            msg_out = [msg_in[i] if i < msg_in_size else 0 for i in range(msg_out_size)]

            quotient, remainder = gf_poly_div(msg_out, gen)

            # Reverse the remainder to match the order used in the original C++ code
            remainder.reverse()

            # Update msg_out with the encoded message
            msg_out = [msg_in[i] if i < msg_in_size else remainder[msg_out_size - 1 - i] for i in range(msg_out_size)]

            return msg_out
        else:
            raise ValueError("The total number of characters - messages + redundant code - exceeds 256")

    except ValueError as ex:
        print("\n\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
        print("Error:", ex)
        print("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")

