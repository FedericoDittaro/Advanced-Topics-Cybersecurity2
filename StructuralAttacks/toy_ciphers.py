s_3b = [0x6, 0x4, 0x5, 0x0, 0x7, 0x1, 0x3, 0x2]
s_4b = [0xe, 0x5, 0xd, 0x1, 0x2, 0xf, 0xc, 0x8, 0x3, 0xa, 0x6, 0xb, 0x4, 0x9, 0x0, 0x7]
p_8b = [0,4,1,5,2,6,3,7]
p_16b = [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]


k_0_3b = 0x4
k_1_3b = 0x7

k_0_4b = 0x3
k_1_4b = 0xa

k_0_8b = [0x8, 0x2]
k_1_8b = [0x1, 0x7]

k_2_8b = [0x3, 0xA]
k_3_8b = [0x6, 0xC]



def encrypt_sub_1r(state):
    state = state ^ k_0_3b
    state = s_3b[state]
    state = state ^ k_1_3b

    return state


def decrypt_sub_1r(state):
    inv_s_box = [s_3b.index(i) for i in range(len(s_3b))]
    state = state ^ k_1_3b
    state = inv_s_box[state]
    state = state ^ k_0_3b

    return state


def encrypt_sub_2r(state):
    state = state ^ k_0_4b
    state = s_4b[state]
    state = state ^ k_1_4b
    state = s_4b[state]

    return state


def decrypt_sub_2r(state):
    inv_s_box = [s_4b.index(i) for i in range(len(s_4b))]
    state = inv_s_box[state]
    state = state ^ k_1_4b
    state = inv_s_box[state]
    state = state ^ k_0_4b

    return state


def encrypt_sub_2r_spn(state): 
    state = [(state >> 4) & 0xF, state & 0xF]

    state = [state[i] ^ k_0_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]
    binary_str = ''.join(f"{x:04b}" for x in state)
    permuted_str = ''.join(binary_str[i] for i in p_8b)
    state = [int(permuted_str[:4], 2), int(permuted_str[4:], 2)]

    state = [state[i] ^ k_1_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]

    return (state[0] << 4) | state[1]


def decrypt_sub_2r_spn(state): 
    inv_s_box = [s_4b.index(i) for i in range(len(s_4b))]
    inv_p_8box = [p_8b.index(i) for i in range(len(p_8b))]

    state = [(state >> 4) & 0xF, state & 0xF]

    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_1_8b[i] for i in range(len(state))]

    binary_str = ''.join(f"{x:04b}" for x in state)
    reversed_str = ''.join(binary_str[inv_p_8box[i]] for i in range(len(inv_p_8box)))
    state = [int(reversed_str[i:i+4], 2) for i in range(0, len(reversed_str), 4)]
    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_0_8b[i] for i in range(len(state))]

    return (state[0] << 4) | state[1]


def encrypt_sub_4r_spn(state):  
    state = [(state >> 4) & 0xF, state & 0xF]

    state = [state[i] ^ k_0_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]
    binary_str = ''.join(f"{x:04b}" for x in state)
    permuted_str = ''.join(binary_str[i] for i in p_8b)
    state = [int(permuted_str[:4], 2), int(permuted_str[4:], 2)]

    state = [state[i] ^ k_1_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]
    binary_str = ''.join(f"{x:04b}" for x in state)
    permuted_str = ''.join(binary_str[i] for i in p_8b)
    state = [int(permuted_str[:4], 2), int(permuted_str[4:], 2)]

    state = [state[i] ^ k_2_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]
    binary_str = ''.join(f"{x:04b}" for x in state)
    permuted_str = ''.join(binary_str[i] for i in p_8b)
    state = [int(permuted_str[:4], 2), int(permuted_str[4:], 2)]

    state = [state[i] ^ k_3_8b[i] for i in range(len(state))]
    state = [s_4b[x] for x in state]

    return (state[0] << 4) | state[1]


def decrypt_sub_4r_spn(state): 
    inv_s_box = [s_4b.index(i) for i in range(len(s_4b))]
    inv_p_8box = [p_8b.index(i) for i in range(len(p_8b))]

    state = [(state >> 4) & 0xF, state & 0xF]

    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_3_8b[i] for i in range(len(state))]

    binary_str = ''.join(f"{x:04b}" for x in state)
    reversed_str = ''.join(binary_str[inv_p_8box[i]] for i in range(len(inv_p_8box)))
    state = [int(reversed_str[i:i+4], 2) for i in range(0, len(reversed_str), 4)]
    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_2_8b[i] for i in range(len(state))]

    binary_str = ''.join(f"{x:04b}" for x in state)
    reversed_str = ''.join(binary_str[inv_p_8box[i]] for i in range(len(inv_p_8box)))
    state = [int(reversed_str[i:i+4], 2) for i in range(0, len(reversed_str), 4)]
    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_1_8b[i] for i in range(len(state))]

    binary_str = ''.join(f"{x:04b}" for x in state)
    reversed_str = ''.join(binary_str[inv_p_8box[i]] for i in range(len(inv_p_8box)))
    state = [int(reversed_str[i:i+4], 2) for i in range(0, len(reversed_str), 4)]
    state = [inv_s_box[x] for x in state]
    state = [state[i] ^ k_0_8b[i] for i in range(len(state))]

    return (state[0] << 4) | state[1]