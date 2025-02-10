import h5py
import numpy as np
from scipy.stats import norm, multivariate_normal
np.set_printoptions(threshold=np.inf)

HW = np.array([0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3,
               3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
               3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2,
               2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5,
               3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5,
               5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3,
               2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
               4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
               3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4,
               4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
               5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5,
               5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8])

AES_SBOX = [
    [0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76],
    [0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0],
    [0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15],
    [0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75],
    [0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84],
    [0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf],
    [0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8],
    [0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2],
    [0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73],
    [0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb],
    [0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79],
    [0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08],
    [0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a],
    [0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e],
    [0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf],
    [0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16]
]

def knownPCor(pred_leak, traces):

    pcor_dist = np.zeros((1, 5))
    # apply a Hamming weight leakage model
    for i in range(5):
        cor = np.corrcoef(pred_leak.T, traces[:, i], rowvar=False)
        pcor_dist[0,i] = cor[0, 1]
    return pcor_dist



if __name__ == '__main__':

    # decide which methods you wish to execute

    np.set_printoptions(precision=4,suppress=True)
    num_points = 700

    # load the dataset
    # power traces and input bytes are stored in a HDF5 file
    # first column corresponds to a single byte, which is the first AES state byte after key XOR: s[0]+k[0] = s[0] + 0
    # second to sixth columns correspond to the actual power trace

    filename = "./ASCAD.h5"

    fhandle = h5py.File(filename)
    profiling_traces = np.array(fhandle['Profiling_traces']['traces'][:, :num_points])
    profiling_labels = np.array(fhandle['Profiling_traces']['labels'], dtype='int')

    attack_traces = np.array(fhandle['Attack_traces']['traces'][:, :num_points])
    attack_labels = np.array(fhandle['Attack_traces']['labels'], dtype='int')

    # For each intermediate calculate noise and signal variance
    elNoise  = np.zeros((1, num_points))
    elSignal = np.zeros((1, num_points))
    snrtrace = np.zeros((1, num_points))

    # Recall that each trace point consists of part that is data dependent (signal)
    # and one part that is data independent (noise)
    # We wrote this as t = L(x) + R


    # We don't know if the signal is related to the HW of the intermediate value s[0]+k[0], but we can test for
    # this via computing the SNR of each trace point under this assumption
    # So we are testing if t = a*HW(x) + r by grouping trace according to HW(x)

    # Expected variance of noise, when traces are grouped by HW
    el = np.zeros((9, num_points))
    for i in range(9):
        ind = np.nonzero(HW[profiling_labels] == i)[0] 
        el[i, :] = np.var(profiling_traces[ind, :], axis=0) 
    elNoise = np.mean(el, axis=0)
    print("[+] ")
    print("The noise in each trace point is:")
    print("%s" % np.array_str(elNoise))
    print("[+] ")
    # Variance of signal, when traces are grouped by HW
    me = np.zeros((9, num_points))
    for h in range(9):
        ind = np.nonzero(HW[profiling_labels] == h)
        me[h, :] = np.mean(profiling_traces[ind, :], axis=1)
    elSignal = np.var(me, axis=0)

    print("[+] ")
    print("The signal in each trace point is:")
    print("%s" % np.array_str(elSignal))
    print("[+] ")

    snrtrace = np.true_divide(elSignal, elNoise)
    print("[+] ")
    print("The SNR in each trace point is:")
    print("%s" % np.array_str(snrtrace))
    print("[+] ")

    # Notice that in order to do the SNR computation we estimated a leakage model for each HW group.
    # This may sound strange, but notice that we assumed t = a*HW(x) + R, so in effect we are trying
    # to find out what the value of the variable a is
    print("[+] ")
    print("The estimated mean for each HW for the last trace point is:")
    print("%s" % np.array_str(me[:,5]))
    print("[+] ")
    # We can test the assumption about the HW being a good model also via doing a correlation analysis with the known key value

    pcor_vals = knownPCor(HW[profiling_labels], profiling_traces)
    print("[+] ")
    print("The Pearson correlation for each trace point is:")
    print("%s" % np.array_str(pcor_vals))

    # Create Gaussian Templates
    print("[+] Creating Gaussian Templates")
    mean_vectors = []
    covariance_matrices = []

    for h in range(9):
        indices = np.nonzero(HW[profiling_labels] == h)
        group_traces = profiling_traces[indices, :].squeeze()

        mean_vector = np.mean(group_traces, axis=0)
        covariance_matrix = np.cov(group_traces.T)

        #determinant = np.linalg.det(covariance_matrix)
        #print(f"Determinant of covariance matrix for HW {h}: {determinant}")
 
        mean_vectors.append(mean_vector)
        covariance_matrices.append(covariance_matrix)

    print("[+] Gaussian Templates Created")

    num_traces = 1000
    attack_traces = attack_traces[:num_traces]
    attack_labels = attack_labels[:num_traces]

    print("[+] Testing Templates with Attack Traces")
    recovered_key = []

    for key_byte_index in range(16):
        print(f"Testing for key byte index: {key_byte_index}")
        predicted_keys = []

        for trace in attack_traces:
            key_likelihoods = []

            for key_guess in range(256): 

                sbox_output = key_guess ^ attack_labels[key_byte_index]
                hamming_weight = HW[sbox_output]

                mean_vector = mean_vectors[hamming_weight]
                covariance_matrix = covariance_matrices[hamming_weight]

                likelihood = multivariate_normal.pdf(trace, mean=mean_vector, cov=covariance_matrix, allow_singular=True)
                key_likelihoods.append(likelihood)

            best_key = np.argmax(key_likelihoods)
            predicted_keys.append(best_key)

        predicted_key_byte = np.argmax(np.bincount(predicted_keys))
        recovered_key.append(predicted_key_byte)

    recovered_key = np.array(recovered_key)
    print("[+] Recovered AES Key:")
    print(recovered_key)