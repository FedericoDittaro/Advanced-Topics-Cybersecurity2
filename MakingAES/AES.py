import math
from pprint import pprint

sbox = [
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16,
]


def GF28_add(x, y):                  # Addition corresponds to bit-wise addition mod 2
    return x^y

def GF28_sub(x, y):                  # Subtraction corresponds to bit-wise addition mod 2
    return x^y

def GF28_multiply(x, y):
    p = 0b100011011             # Using the AES irreducible polynomial x^8+x^4+x^3+x+1
    m = 0                       # Performing school book multiplication (bit-wise), m holds product
    for i in range(8):
        m = m << 1              # Left shift intermediate sum each bit
        if m & 0b100000000:     # If larger than 255 then reduce
            m = m ^ p
        if y & 0b010000000:     # If multiplier bit is set then add y
            m = m ^ x
        y = y << 1              # Left shift multiplier to check next bit in next iteration
    return m

def xTimes(x):               # Your turn: implement the special case where we multiply by "x" (i.e. 2) without doing an explicit finite field multiplication
    x = x << 1
    if x & 0b100000000:
        x = x ^ 0b100011011

    return x


def GF28_inv(x):                    # Your turn: implement finite field inversion
    if x == 0:
        return 0

    def gf_pow(base, exp):
        result = 1
        while exp > 0:
            if exp & 1:
                result = GF28_multiply(result, base)
            base = GF28_multiply(base, base)
            exp >>= 1
        return result

    return gf_pow(x, 254)


def Sbox_generation(byte):
    c = 0b01100011
    inverse = GF28_inv(byte)
    new_byte = 0
    for i in range(8):
        bit = (inverse >> i) & 1 
        transformed_bit = bit ^ ((inverse >> ((i+4) % 8)) & 1) \
                               ^ ((inverse >> ((i+5) % 8)) & 1) \
                               ^ ((inverse >> ((i+6) % 8)) & 1) \
                               ^ ((inverse >> ((i+7) % 8)) & 1) \
                               ^ ((c >> i) & 1) 
        new_byte |= (transformed_bit << i)
    
    return new_byte

def SubBytes(state):
    return [[Sbox_generation(byte) for byte in row] for row in state]


def pprint_state(state):
    pprint([["{0:0{1}x}".format(i, 2).upper() for i in row] for row in state])


#def SubBytes(state):
#    return [[sbox[i] for i in row] for row in state]

def ShiftRows(state):
    return [
        state[0],
        state[1][1:] + state[1][:1],
        state[2][2:] + state[2][:2],
        state[3][3:] + state[3][:3],
    ]

def SubBytesAndShiftRows(state):
    return [
        [Sbox_generation(byte) for byte in state[0]],
        [Sbox_generation(byte) for byte in (state[1][1:] + state[1][:1])],  
        [Sbox_generation(byte) for byte in (state[2][2:] + state[2][:2])],  
        [Sbox_generation(byte) for byte in (state[3][3:] + state[3][:3])],  
    ]


def AddRoundKey(state, rkey):
    pairs = [zip(a, b) for a, b in zip(state, rkey)]
    return [[a ^ b for a, b in row] for row in pairs]

def MixColumns(state):                                  # Your turn to implement MixColumn
    def mix_single_column(column):
        return [
            xTimes(column[0]) ^ GF28_multiply(column[1], 3) ^ column[2] ^ column[3],
            column[0] ^ xTimes(column[1]) ^ GF28_multiply(column[2], 3) ^ column[3],
            column[0] ^ column[1] ^ xTimes(column[2]) ^ GF28_multiply(column[3], 3),
            GF28_multiply(column[0], 3) ^ column[1] ^ column[2] ^ xTimes(column[3])
        ]

    transposed_state = [list(col) for col in zip(*state)]
    mixed_columns = [mix_single_column(col) for col in transposed_state]
    return [list(row) for row in zip(*mixed_columns)]


if __name__ == "__main__":

    # Consider the state to be two-dimensional array of bytes
    # The order is defined in FIPS 197
    AES_state = [
        [0x32, 0x88, 0x31, 0xe0],
        [0x43, 0x5a, 0x31, 0x37],
        [0xf6, 0x30, 0x98, 0x07],
        [0xa8, 0x8d, 0xa2, 0x34],
    ]

    k1 = [
        [0x2b, 0x28, 0xab, 0x09],
        [0x7e, 0xae, 0xf7, 0xcf],
        [0x15, 0xd2, 0x15, 0x4f],
        [0x16, 0xa6, 0x88, 0x3c],
    ]
    print('AES State beginning of round')
    pprint_state(AES_state)

    print('AES State after AddRoundKey')
    state = AddRoundKey(AES_state,k1)
    pprint_state(state)

    #print('AES State after SubBytes')
    #state = SubBytes(state)
    #pprint_state(state)

    #print('AES State after ShiftRows')
    #state = ShiftRows(state)
    #pprint_state(state)

    print('AES State using the combination of SubBytes and ShiftRows')
    state = SubBytesAndShiftRows(state)
    pprint_state(state)

    print('AES State after MixColumns')
    state = MixColumns(state)
    pprint_state(state)