import random
from typing import List

from gf import *
from help_gf import *
from rs_decoding import *
from rs_encoding import *


D = 9
N = 31
k = N - D + 1

message = b"Hi my freinds! How are you?"
msglen = len(message)

print("Messege: \n")

print(message)

print("Original:  ")
original = bytearray()



for i in range(msglen):
    original.append(message[i])
    print(original[i], end=" ")







encoded = rs_encode_msg(original, N - k)

print("\nEncoded:   ")
for i in range(len(encoded)):
    print(encoded[i], end=" ")



encoded_str = [x for x in encoded]
erroneous = list(encoded_str)


print("\nErroneous: ")


random.seed()
for i in range(len(encoded)):
    if i < (N - k) // 2:
        h = random.randint(0, msglen - 1)
        r = random.randint(0, 255)
        erroneous[h] = r


for i in range(len(erroneous)):
    print(erroneous[i], end=" ")




decoded = rs_decode_msg(erroneous,N-k)             #N-k
print("\n")

print("decoded: ")
for i in range(len(decoded)):
    print(decoded[i], end=" ")


print("\n")
for i in range(min(msglen, len(decoded))):
    char_value = decoded[i]
    char_symbol = chr(char_value)
    print(char_symbol,end="")







