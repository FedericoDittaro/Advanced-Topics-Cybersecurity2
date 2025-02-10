from pprint import pprint
import AES

def GF28_multiply(x, y):
    p = 0b100011011             # Using the AES irreducible polynomial x^8+x^4+x^3+x+1
    m = 0                       # Performing school book multiplication (bit-wise), m holds product
    for i in range(8):
        m = m << 1              # Left shift intermediate sum each bit
        if m & 0b100000000:     # If larger than 255 then reduce
            pprint("At step ", i, "reduction performed")
            m = m ^ p
        if y & 0b010000000:     # If multiplier bit is set then add y
            #pprint("At step ", i, "m XOR x")
            m = m ^ x
            #pprint("Value of m ", m)
            
        y = y << 1              # Left shift multiplier to check next bit in next iteration
    return m


def xTimes(n, x):
    def multiply_by_2(val):
        val = val << 1
        if val & 0b100000000:
            val ^= 0b100011011
        return val

    m = 0
    p = n

    while x > 0:
        if x & 1:
            m = m ^ p
        p = multiply_by_2(p)
        x >>= 1

    return m


if __name__ == "__main__":

   ris = 0b00010110 ^ 0b00001011
   m = 11 
   pprint(ris)