from dnds import *
from helperfuncs import *
import string
import random
from matplotlib import pyplot as plt

if __name__ == "__main__":
    # generate 10k long genome string
    num_runs = 100
    str_len = 10000
    mut_rate = 0.01
    seed = 0
    k = 7
    scale_factor = 0.3
    random.seed(seed)

    for scale_factor in [0.1, 0.2, 0.3, 0.5, 0.75, 1.0]:
        correct_dnds_values = []
        estimated_dnds_values_whole_seq = []
        estimated_dnds_values_fmh = []
        for i in range(num_runs):
            s1 = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))
            s2 = mutate_string(s1, mut_rate)

            approx_dnds_using_whole_seq = determine_approx_dnds_no_fmh(s1, s2, k)
            correct_dnds = determine_correct_dnds(s1, s2)
            approx_dnds_using_fmh = determine_approx_dnds_using_fmh(s1, s2, k, scale_factor)

            correct_dnds_values.append(correct_dnds)
            estimated_dnds_values_whole_seq.append(approx_dnds_using_whole_seq)
            estimated_dnds_values_fmh.append(approx_dnds_using_fmh)

        plt.clf()
        plt.scatter(correct_dnds_values, estimated_dnds_values_whole_seq, label='Using whole seq', alpha=0.9)
        plt.scatter(correct_dnds_values, estimated_dnds_values_fmh, label='Using FMH with scalef = ' + str(scale_factor), alpha=0.5)
        tmp = [min(correct_dnds_values), max(correct_dnds_values)]
        plt.plot(tmp, tmp, linestyle='--')
        plt.ylim(0,20)
        plt.xlabel('Correct dnds value')
        plt.ylabel('Estimated dnds value')
        plt.legend()
        plt.savefig('approx_vs_real_dnds_scale_' + str(scale_factor) + '.pdf')
