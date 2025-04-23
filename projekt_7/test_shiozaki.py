from GF_arithmetic import *
from RS_encoder import rs_encode_msg
from RS_decoder_error_locator import rs_shiozaki_gao_decode
import random

def main():
    prim = 0x11d
    n = 20
    k = 11
    nsym = n - k
    message = "hello world"  # 11 znaków

    init_tables(prim)

    msg_bytes = [ord(c) for c in message]
    print("Wiadomość:", msg_bytes)

    codeword = rs_encode_msg(msg_bytes, nsym)
    print("Zakodowana wiadomość:", codeword)

    corrupted = codeword[:]
    for i in [0, 1, 3]:
        corrupted[i] ^= random.randint(1, 255)
    print("Uszkodzona wiadomość:", corrupted)

    decoded_msg, ecc = rs_shiozaki_gao_decode(corrupted, nsym)
    print("Zdekodowana wiadomość:", decoded_msg)
    print("Zdekodowana jako tekst:", ''.join(chr(x) for x in decoded_msg))


if __name__ == "__main__":
    main()
