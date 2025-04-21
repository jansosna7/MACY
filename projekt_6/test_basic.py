#https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
from GF_arithmetic import *
from RS_encoder import *
from RS_decoder_error_locator import *

import random
# Configuration of the parameters and input message
prim = 0x11d
n = 20 # set the size you want, it must be > k, the remaining n-k symbols will be the ECC code (more is better)
k = 11 # k = len(message)
message = "hello world" # input message

# Initializing the log/antilog tables
init_tables(prim)

# Encoding the input message
mesecc = rs_encode_msg([ord(x) for x in message], n-k)
print("Original: %s" % mesecc)

# Tampering 6 characters of the message (over 9 ecc symbols, so we are above the Singleton Bound)

for i in range(6):
    mesecc[random.randint(0,k-1)] = random.randint(0,255)
    
print("Corrupted: %s" % mesecc)

# Decoding/repairing the corrupted message, by providing the locations of a few erasures, we get below the Singleton Bound
# Remember that the Singleton Bound is: 2*e+v <= (n-k)
corrected_message, corrected_ecc = rs_correct_msg(mesecc, n-k, erase_pos=[0, 1, 2])
print("Repaired: %s" % (corrected_message+corrected_ecc))
print(''.join([chr(x) for x in corrected_message]))
