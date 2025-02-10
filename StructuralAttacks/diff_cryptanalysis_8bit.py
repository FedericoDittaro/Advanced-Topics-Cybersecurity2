import random
import toy_ciphers
from pprint import pprint

# Compute the difference distribution table for a specific sbox
def difference_distribution_table(sbox):

    d_dist_table= [[0 for row in range(len(sbox))] for col in range(len(sbox))]
    for x in range(len(sbox)):
        for x_d in range(len(sbox)):
            d_in = x ^ x_d

            y = sbox[x]
            y_d = sbox[x_d]
            d_out = y ^ y_d

            d_dist_table[d_in][d_out] = d_dist_table[d_in][d_out] +1

    return d_dist_table

def get_good_pairs(d_in, d_out, sbox):    # d_in: input difference      d_out: output_difference

    pairs= []

    for x in range(len(sbox)):
        x_d = x ^ d_in
        y = sbox[x]
        y_d = sbox[x_d]
        d_out_test = y ^ y_d
        if d_out_test == d_out:
            pairs.append((x,y))

    return pairs

def inverse_sbox_pbox (m):
    m = [(m >> 4) & 0xF, m & 0xF]
    m = [inv_s_box[x] for x in m]
    binary_str_c = ''.join(f"{x:04b}" for x in m)
    reversed_str_c = ''.join(binary_str_c[inv_p_box[i]] for i in range(len(inv_p_box)))
    m = [int(reversed_str_c[i:i+4], 2) for i in range(0, len(reversed_str_c), 4)]
    return m



if __name__=="__main__":

    # demo for the principle of differential cryptanalysis
    # the demo is for the 8 bit, two round cipher, features two independent subkeys
    # grab the corresponding sbox from the toy_ciphers
    s_box = toy_ciphers.s_4b
    p_box = toy_ciphers.p_8b
    inv_s_box = [s_box.index(i) for i in range(len(s_box))]
    inv_p_box = [p_box.index(i) for i in range(len(p_box))]

    diff_dist_table = difference_distribution_table(s_box)
    print('Difference table:\n')
    pprint(diff_dist_table)
    print('')
    # find a good difference pair via a visual inspection of the table,
    # one such choice is used in the next line
    good_difference = (6,2)
    good_pairs = get_good_pairs(good_difference[0],good_difference[1], s_box)
    print('For the input-output difference',good_difference,'the good (u,v) pairs are:')
    pprint(good_pairs)
    print('')
    # pick a message (at random) and then find the corresponding
    # second message via a good differential, repeat until we find
    # a good output differential
    good_input_right=[]
    good_input_left=[]
    for m_right in range(16):
        m_d_right = m_right^6
        m_left = m_right << 4
        m_d_left = m_left^(6 << 4)

        c_right =  toy_ciphers.encrypt_sub_2r_spn(m_right)
        c_right = inverse_sbox_pbox(c_right)
        c_d_right = toy_ciphers.encrypt_sub_2r_spn(m_d_right)
        c_d_right = inverse_sbox_pbox(c_d_right)
        c_left = toy_ciphers.encrypt_sub_2r_spn(m_left)
        c_left = inverse_sbox_pbox(c_left)
        c_d_left = toy_ciphers.decrypt_sub_2r_spn(m_d_left)
        c_d_left = inverse_sbox_pbox(c_d_left)

        if c_right[1]^c_d_right[1] ==2:
            good_input_right.append(((m_right,m_d_right),(c_right[1],c_d_right[1])))

        if c_left[0]^c_d_left[0] ==2:
            good_input_left.append(((m_left >> 4,m_d_left >> 4),(c_left[0],c_d_left[0])))
    print('For the above input-output difference, the good input-output pairs ((m, m+delta_m),(c, c+delta_c)) are:')
    pprint(good_input_right)
    pprint(good_input_left)
    print('')


    # pick a good pair (or try as many as there are according to the selected differential)
    # and work out the keys

    # we know that m + k0 must be one of the inputs to the sbox according to the differential
    # (and there are only a few according to our characteristic that are possible)
    # we also know that c+k1 must be one of the outputs of the sbox accord. to the differential
    k0_list_right: list[int] = [0]*len(good_pairs)
    k0_list_left: list[int] = [0]*len(good_pairs)
    k1_list_right: list[int] = [0]*len(good_pairs)
    k1_list_left: list[int] = [0]*len(good_pairs)


    for i in range(len(good_pairs)):
        k0_list_right[i]= good_pairs[i][0] ^ good_input_right[0][0][0]
        k0_list_left[i]= good_pairs[i][0] ^ good_input_left[0][0][0]
        k1_list_right[i]= good_pairs[i][1] ^ good_input_right[0][1][0]
        k1_list_left[i]= good_pairs[i][1] ^ good_input_left[0][1][0]


    print('Candidates for k0 left:', k0_list_left)
    print('Candidates for k0 right:', k0_list_right)
    print('Candidates for k1 left:', k1_list_left)
    print('Candidates for k1 right:', k1_list_right)
    message = random.randint(0, 15)
    for i in range(len(k0_list_left)):
        for j in range(len(k0_list_right)):

            cipher = toy_ciphers.encrypt_sub_2r_spn(message)
            challenge = [(message >> 4) & 0xF, message & 0xF]

            challenge[0] = challenge[0] ^ k0_list_left[i]
            challenge[1] = challenge[1] ^ k0_list_right[j]
            challenge = [s_box[x] for x in challenge]
            binary_str = ''.join(f"{x:04b}" for x in challenge)
            permuted_str = ''.join(binary_str[i] for i in p_box)
            challenge = [int(permuted_str[:4], 2), int(permuted_str[4:], 2)]

            challenge[0] = challenge[0] ^ k1_list_left[i]
            challenge[1] = challenge[1] ^ k1_list_right[j]
            challenge = [s_box[x] for x in challenge]
            challenge = (challenge[0] << 4) | challenge[1]


            if cipher == challenge:
                print('Message: ', message)
                print('Keys:', '[', k0_list_left[i], ',', k0_list_right[j], ']', ' [', k1_list_left[i], ',', k1_list_right[j], ']')