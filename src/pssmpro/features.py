import re
import os
import numpy as np
import itertools

# Global Variable, all available encoders
all_encoders = ['aac_pssm', 'aadp_pssm', 'aatp', 'ab_pssm', 'd_fpssm', 'dp_pssm', 'dpc_pssm', 'edp', 'eedp',
                'k_separated_bigrams_pssm', 'medp', 'pse_pssm', 'pssm_ac', 'pssm_cc',
                'pssm_composition', 'rpm_pssm', 'rpssm', 's_fpssm', 'smoothed_pssm', 'tpc', 'tri_gram_pssm']


# pssm input file reading functions
def read_pssm_matrix(input_matrix):
    """ Reading the PSSM input file into a numpy array"""
    PSSM = []
    p = re.compile(r'-*[0-9]+')
    stream = open(input_matrix)
    for line, string in enumerate(stream.readlines()):
        if line > 2:
            str_vec = []
            overall_vec = string.split()
            if len(overall_vec) == 0:
                break
            str_vec.extend(overall_vec[1])

            if len(overall_vec) < 44:
                for cur_str in overall_vec[2:]:
                    str_vec.extend(p.findall(cur_str))
                    if len(str_vec) >= 21:
                        if len(str_vec) > 21:
                            raise ValueError("Wrong PSSM format")
                        break
                print("Done")
            else:
                str_vec = string.split()[1:42]
            if len(str_vec) == 0:
                break
            PSSM.append(str_vec)
    PSSM = np.array(PSSM)
    return PSSM


# pssm matrix parsing functions
def average(matrix_sum, seq_len):
    """average the rows of a pssm"""
    matrix_array = np.divide(matrix_sum, seq_len)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1],))]
    return matrix_average


def normalize_pssm(pssm):
    """ Normalizing a PSSM """
    pssm = pssm[:, 1:21]
    pssm = pssm.astype(float)

    seq_cn = np.shape(pssm)[0]
    PSSM_norm = [[0.0] * 20] * seq_cn
    PSSM_norm = np.array(PSSM_norm)
    mean_matrix = np.mean(pssm, axis=1)
    std_matrix = np.std(pssm, axis=1)

    for i in range(seq_cn):
        for j in range(20):
            if std_matrix[i] == 0.0:
                PSSM_norm[i][j] = pssm[i][j] - mean_matrix[i]
            else:
                PSSM_norm[i][j] = (pssm[i][j] - mean_matrix[i]) / std_matrix[i]
    return PSSM_norm


def window(pssm, w_smooth, w_slide):
    """Smooth PSSM creation"""
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    w_smooth = int(w_smooth)
    w_slide = int(w_slide)

    pssm = pssm[:, 1:21]
    pssm = pssm.astype(float)

    seq_cn = np.shape(pssm)[0]

    # original PSSM
    PSSM_smooth = np.array([[0.0] * 20] * seq_cn)
    PSSM_orig = np.array(pssm)

    # section for PSSM_smooth features
    PSSM_smooth_full = pssm_smooth(PSSM_orig, PSSM_smooth, w_smooth, seq_cn)

    PSSM_smooth_final = [[0.0] * 20] * w_slide
    PSSM_smooth_final = np.array(PSSM_smooth_final)
    for i in range(w_slide):
        PSSM_smooth_final[i] = PSSM_smooth_full[i]
    matrix_final = average(PSSM_smooth_final, 1.0)
    return matrix_final


def pssm_smooth(pssm_original, pssm_smooth_new, w_smooth, seq_len):
    """Smooth PSSM creation helper"""
    for i in range(seq_len):
        if i < (w_smooth - 1) / 2:
            for j in range(i + (w_smooth - 1) // 2 + 1):
                pssm_smooth_new[i] += pssm_original[j]
        elif i >= (seq_len - (w_smooth - 1) // 2):
            for j in range(i - (w_smooth - 1) // 2, seq_len):
                pssm_smooth_new[i] += pssm_original[j]
        else:
            for j in range(i - (w_smooth - 1) // 2, i + (w_smooth - 1) // 2 + 1):
                pssm_smooth_new[i] += pssm_original[j]
    return pssm_smooth_new


def handle_rows(pssm, switch, count):
    """
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.

    Convert the PSSM matrix into sum of their rows either by taking the amino acids or not taking them into account
    """
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    matrix_final = [[0.0] * 20] * int(count / 20)
    matrix_final = np.array(matrix_final)
    seq_cn = 0

    PSSM_shape = np.shape(pssm)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec = pssm[i]
        str_vec_positive = list(map(int, str_vec[1:21]))
        str_vec_positive = np.array(str_vec_positive)
        if switch == 1:
            str_vec_positive[str_vec_positive < 0] = 0
        elif switch == 2:
            str_vec_positive[str_vec_positive < 0] = 0
            str_vec_positive[str_vec_positive > 7] = 0
        if count == 20:
            matrix_final[0] = list(map(sum, zip(str_vec_positive, matrix_final[0])))
        elif count == 400:
            matrix_final[Amino_vec.index(str_vec[0])] = list(map(sum, zip(str_vec_positive,
                                                                          matrix_final[Amino_vec.index(str_vec[0])])))

    return matrix_final


def pre_handle_columns(pssm, step, part, idx):
    """
    if step=k, we calculate the relation between one residue and the kth residue afterward.

    if part=0, we calculate the left part of PSSM.
    if part=1, we calculate the right part of PSSM.

    if idx=0, we product the residue-pair.
    if idx=1, we minus the residue-pair.

    if key=1, we divide each element by the sum of elements in its column.
    if key=0, we don't perform the above process.

    Preprocessing the PSSM column-wise
    """

    if part == 0:
        pssm = pssm[:, 1:21]
    elif part == 1:
        pssm = pssm[:, 21:]
    pssm = pssm.astype(float)
    matrix_final = [[0.0] * 20] * 20
    matrix_final = np.array(matrix_final)
    seq_cn = np.shape(pssm)[0]

    if idx == 0:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - step):
                    matrix_final[i][j] += (pssm[k][i] * pssm[k + step][j])

    elif idx == 1:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - step):
                    matrix_final[i][j] += ((pssm[k][i] - pssm[k + step][j]) * (pssm[k][i] - pssm[k + step][j]) / 4.0)

    return matrix_final


def handle_tri_columns(pssm):
    """
    Preprocessing three columns of a PSSM
    """
    matrix_final = [[[0.0] * 20] * 20] * 20
    matrix_final = np.array(matrix_final)
    pssm = pssm[:, 21:]
    pssm = pssm.astype(float)
    pssm = np.array(pssm)

    seq_cn = np.shape(pssm)[0]
    for m in range(20):
        for n in range(20):
            for r in range(20):
                for i in range(seq_cn - 2):
                    matrix_final[m][n][r] += (pssm[i, m] * pssm[i + 1, n] * pssm[i + 2, r])

    matrix_final = np.divide(matrix_final, 1000000.0)
    matrix_final_shape = np.shape(matrix_final)

    matrix_result = [
        (np.reshape(matrix_final, (matrix_final_shape[0] * matrix_final_shape[1] * matrix_final_shape[2],)))]
    return matrix_result


def handle_mixed(pssm, alpha):
    """
    preprocessing PSSM mixed
    """
    row1 = [0.0] * 20
    row2 = [0.0] * 20

    matrix_final = [[0.0] * 40] * 1
    row1 = np.array(row1)
    row2 = np.array(row2)
    matrix_final = np.array(matrix_final)

    pssm_norm = normalize_pssm(pssm)
    seq_cn = np.shape(pssm)[0]
    for i in range(seq_cn):
        row1 = list(map(sum, zip(row1, pssm_norm[i])))

    row1 = np.divide(row1, seq_cn)

    for j in range(20):
        for i in range(seq_cn - alpha):
            row2[j] += (pssm_norm[i][j] - pssm_norm[i + alpha][j]) * (pssm_norm[i][j] - pssm_norm[i + alpha][j])

    row2 = np.divide(row2, seq_cn - alpha)

    row = np.hstack((row1, row2))
    matrix_final[0] = row
    return matrix_final


def handle_mixed2(pssm, alpha):
    """
    preprocessing PSSM mixed2
    :param pssm:
    :param alpha:
    :return:
    """
    row1 = [0.0] * 40
    row2 = [[0.0] * (2 * alpha)] * 20
    matrix_final = [[0.0] * (40 + 40 * alpha)] * 1

    row1 = np.array(row1)
    row2 = np.array(row2)
    matrix_final = np.array(matrix_final)

    pssm_norm = normalize_pssm(pssm)
    seq_cn = np.shape(pssm)[0]
    for j in range(20):
        positive_count_1 = 0
        negative_count_1 = 0
        for i in range(seq_cn):
            if pssm_norm[i][j] >= 0:
                positive_count_1 += 1
                row1[2 * j] += pssm_norm[i][j]
            elif pssm_norm[i][j] < 0:

                negative_count_1 += 1
                row1[2 * j + 1] += pssm_norm[i][j]

        row1[2 * j] = row1[2 * j] / positive_count_1
        row1[2 * j + 1] = row1[2 * j + 1] / negative_count_1

    for j in range(20):
        for alpha in range(1, alpha + 1):
            positive_count_2 = 0
            negative_count_2 = 0
            for i in range(seq_cn - alpha):
                if (pssm_norm[i][j] - pssm_norm[i + alpha][j]) >= 0:
                    positive_count_2 += 1
                    row2[j][2 * alpha - 2] += (pssm_norm[i][j] - pssm_norm[i + alpha][j]) * (
                                pssm_norm[i][j] - pssm_norm[i + alpha][j])
                elif (pssm_norm[i][j] - pssm_norm[i + alpha][j]) < 0:
                    negative_count_2 += 1
                    row2[j][2 * alpha - 1] += (pssm_norm[i][j] - pssm_norm[i + alpha][j]) * (
                                pssm_norm[i][j] - pssm_norm[i + alpha][j])
            row2[j][2 * alpha - 2] = row2[j][2 * alpha - 2] / positive_count_2
            row2[j][2 * alpha - 1] = row2[j][2 * alpha - 1] / negative_count_2

    row2 = average(row2, 1.0)
    row = np.hstack((row1, row2[0]))
    matrix_final[0] = row
    return matrix_final


def handle_mixed3(pssm):
    """
    preprocessing PSSM mixed3
    """
    row = [[0.0] * 110] * 1
#    row1 = [0.0] * 100
    row2 = [0.0] * 10
    row = np.array(row)
    row2 = np.array(row2)

    seq_cn = np.shape(pssm)[0]
    r_pssm = [[0.0] * 10] * seq_cn
    r_pssm = np.array(r_pssm)

    pssm = pssm[:, 1:21]
    pssm = pssm.astype(float)
    pssm = np.array(pssm)

    r_pssm[:, 0] = np.divide(list(map(sum, zip(pssm[:, 13], pssm[:, 17], pssm[:, 18]))), 3.0)
    r_pssm[:, 1] = np.divide(list(map(sum, zip(pssm[:, 10], pssm[:, 12]))), 2.0)
    r_pssm[:, 2] = np.divide(list(map(sum, zip(pssm[:, 9], pssm[:, 19]))), 2.0)
    r_pssm[:, 3] = np.divide(list(map(sum, zip(pssm[:, 0], pssm[:, 15], pssm[:, 16]))), 3.0)
    r_pssm[:, 4] = np.divide(list(map(sum, zip(pssm[:, 2], pssm[:, 8]))), 2.0)
    r_pssm[:, 5] = np.divide(list(map(sum, zip(pssm[:, 5], pssm[:, 6], pssm[:, 3]))), 3.0)
    r_pssm[:, 6] = np.divide(list(map(sum, zip(pssm[:, 1], pssm[:, 11]))), 2.0)
    r_pssm[:, 7] = pssm[:, 4]
    r_pssm[:, 8] = pssm[:, 7]
    r_pssm[:, 9] = pssm[:, 14]

    mean_matrix = np.mean(r_pssm, axis=0)

    for j in range(10):
        for i in range(seq_cn):
            row2[j] += (r_pssm[i][j] - mean_matrix[j]) * (r_pssm[i][j] - mean_matrix[j])

    row2 = np.divide(row2, seq_cn)
    matrix_final = [[0.0] * 10] * 10
    matrix_final = np.array(matrix_final)
    for i in range(10):
        for j in range(10):
            for k in range(seq_cn - 1):
                matrix_final[i][j] += ((r_pssm[k][i] - r_pssm[k + 1][j]) * (r_pssm[k][i] - r_pssm[k + 1][j]) / 2.0)

    row1 = average(matrix_final, seq_cn - 1)[0]
    row[0] = np.hstack((row1, row2))
    return row


def correlation(pssm, idx, group):
    """
    PSSM correlation calculator
    :param pssm:
    :param idx:
    :param group:
    :return:
    """
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    # GROUP=10
    pssm = pssm[:, 1:21]
    pssm = pssm.astype(float)

    # section for PSSM_AC features
    seq_cn = np.shape(pssm)[0]
    g = group
    seq_len = seq_cn
    if idx == 0:
        matrix_final = pssm_ac_cal(pssm, g, seq_len)
    elif idx == 1:
        matrix_final = pssm_cc_cal(pssm, g, seq_len)
    else:
        raise ValueError("Wrong idx assigned")

    matrix_final = average(matrix_final, seq_len)
    return matrix_final


def pssm_ac_cal(pssm, g, seq_len):
    """
    PSSM AC helper
    :param pssm:
    :param g:
    :param seq_len:
    :return:
    """
    pssm_ac_matrix = np.array([[0.0] * 20] * g)

    for pg in range(g):
        l_g = seq_len - pg - 1
        for pj in range(20):
            sum_jl = 0.0
            for i in range(seq_len):
                sum_jl += pssm[i][pj]
            sum_jl /= seq_len

            pssm_acjg = 0.0
            for i in range(l_g):
                pssm_acjg += (pssm[i][pj] - sum_jl) * (pssm[i + pg + 1][pj] - sum_jl)
            pssm_acjg /= l_g
            pssm_ac_matrix[pg][pj] = pssm_acjg
    return pssm_ac_matrix


def pssm_cc_cal(pssm, g, seq_len):
    """
    PSSM CC helper
    :param pssm:
    :param g:
    :param seq_len:
    :return:
    """
    PSSM_CC = np.array([[0.0] * 380] * g)
    for pg in range(g):
        l_g = seq_len - pg - 1
        for pj_1 in range(20):
            sum_jl_1 = 0.0
            for i in range(seq_len):
                sum_jl_1 += pssm[i][pj_1]
            sum_jl_1 /= seq_len

            for pj_2 in range(20):
                if pj_2 != pj_1:
                    sum_jl_2 = 0.0
                    for i in range(seq_len):
                        sum_jl_2 += pssm[i][pj_2]
                    sum_jl_2 /= seq_len
                    pssm_acjg = 0.0
                    for i in range(l_g):
                        pssm_acjg += (pssm[i][pj_1] - sum_jl_1) * (pssm[i + pg + 1][pj_2] - sum_jl_2)
                    pssm_acjg /= l_g
                    if pj_1 < pj_2:
                        PSSM_CC[pg][19 * pj_1 + (pj_2 - 1)] = pssm_acjg
                    else:
                        PSSM_CC[pg][19 * pj_1 + pj_2] = pssm_acjg

    return PSSM_CC


# Actual functions of encoders
# row transformation based feature encoding algorithms
def aac_pssm(input_matrix):
    """
    AAC PSSM feature encoder
    :param input_matrix:
    :return:
    """
    SWITCH = 0
    COUNT = 20
    seq_cn = float(np.shape(input_matrix)[0])
    aac_pssm_matrix = handle_rows(input_matrix, SWITCH, COUNT)
    aac_pssm_matrix = np.array(aac_pssm_matrix)
    aac_pssm_vector = average(aac_pssm_matrix, seq_cn)
    return aac_pssm_vector[0]


def d_fpssm(input_matrix):
    """
    D FPSSM feature encoder
    :param input_matrix:
    :return:
    """
    SWITCH = 1
    COUNT = 20
    d_fpssm_matrix = handle_rows(input_matrix, SWITCH, COUNT)
    seqLen = float(np.shape(input_matrix)[0])
    seq_cn = 1.0
    element_max = np.amax(d_fpssm_matrix[0], axis=0)
    element_min = np.amin(d_fpssm_matrix[0], axis=0)
    d_fpssm_vector = average(d_fpssm_matrix, seq_cn)
    d_fpssm_vector = np.add(d_fpssm_vector, -element_min)
    d_fpssm_vector = np.divide(d_fpssm_vector, element_max)
    d_fpssm_vector = np.divide(d_fpssm_vector, seqLen)
    return d_fpssm_vector[0]


def smoothed_pssm(input_matrix, smooth=7, slide=50):
    """
    Smoothed PSSM feature encoder
    :param input_matrix:
    :param smooth:
    :param slide:
    :return:
    """
    smoothed_pssm_matrix = window(input_matrix, smooth, slide)
    return smoothed_pssm_matrix[0]


def ab_pssm(input_matrix):
    """
    AB PSSM feature encoder
    :param input_matrix:
    :return:
    """
    seq_cn = np.shape(input_matrix)[0]
    BLOCK = int(seq_cn/20)

    matrix_final = []
    for i in range(19):
        tmp = input_matrix[i*BLOCK:(i+1)*BLOCK]
        matrix_final.append(aac_pssm(tmp))
    tmp = input_matrix[19*BLOCK:]
    matrix_final.append(aac_pssm(tmp))
    ab_pssm_matrix = average(matrix_final, 1.0)
    return ab_pssm_matrix[0]


def pssm_composition(input_matrix):
    """
    PSSM composition feature encoder
    :param input_matrix:
    :return:
    """
    SWITCH = 0
    COUNT = 400
    seq_cn = float(np.shape(input_matrix)[0])
    pssm_composition_matrix = handle_rows(input_matrix, SWITCH, COUNT)
    pssm_composition_vector = average(pssm_composition_matrix, seq_cn)
    return pssm_composition_vector[0]


def rpm_pssm(input_matrix):
    """
    RPM PSSM feature encoder
    :param input_matrix:
    :return:
    """
    SWITCH = 1
    COUNT = 400
    seq_cn = float(np.shape(input_matrix)[0])
    rpm_pssm_matrix = handle_rows(input_matrix, SWITCH, COUNT)
    rpm_pssm_vector = average(rpm_pssm_matrix, seq_cn)
    return rpm_pssm_vector[0]


def s_fpssm(input_matrix):
    """
    S FPSSM feature encoder
    :param input_matrix:
    :return:
    """
    SWITCH = 2
    COUNT = 400
    s_fpssm_matrix = handle_rows(input_matrix, SWITCH, COUNT)
    s_fpssm_matrix = np.array(s_fpssm_matrix)
    s_fpssm_matrix_shape = np.shape(s_fpssm_matrix)
    matrix_average = [(np.reshape(s_fpssm_matrix, (s_fpssm_matrix_shape[0] * s_fpssm_matrix_shape[1], )))]
    return matrix_average[0]


# column transformation based feature encoders
def dpc_pssm(input_matrix):
    """
    DPC PSSM feature encoder
    :param input_matrix:
    :return:
    """
    PART = 0
    STEP = 1
    ID = 0
    matrix_final = pre_handle_columns(input_matrix, STEP, PART, ID)
    seq_cn = float(np.shape(input_matrix)[0])
    dpc_pssm_vector = average(matrix_final, seq_cn-STEP)
    return dpc_pssm_vector[0]


def k_separated_bigrams_pssm(input_matrix, step=1):
    """
    k separated bigrams feature encoder
    :param input_matrix:
    :param step:
    :return:
    """
    PART = 1
    ID = 0
    matrix_final = pre_handle_columns(input_matrix, step, PART, ID)
    k_separated_bigrams_pssm_vector = average(matrix_final, 10000.0)
    return k_separated_bigrams_pssm_vector[0]


def tri_gram_pssm(input_matrix):
    """
    Tri gram PSSM feature encoder
    :param input_matrix:
    :return:
    """
    tri_gram_pssm_matrix = handle_tri_columns(input_matrix)
    return tri_gram_pssm_matrix[0]


def eedp(input_matrix):
    """
    EEDP feature encoder
    :param input_matrix:
    :return:
    """
    STEP = 2
    PART = 0
    ID = 1
    seq_cn = float(np.shape(input_matrix)[0])
    matrix_final = pre_handle_columns(input_matrix, STEP, PART, ID)
    eedp_vector = average(matrix_final, seq_cn-STEP)
    return eedp_vector[0]


def tpc(input_matrix):
    """
    TPC feature encoder
    :param input_matrix:
    :return:
    """
    PART = 0
    STEP = 1
    ID = 0
    matrix_final = pre_handle_columns(input_matrix, STEP, PART, ID)
    matrix_tmp = [0.0] * 20
    matrix_tmp = np.array(matrix_tmp)
    for i in range(20):
        matrix_tmp = list(map(sum, zip(matrix_final[i], matrix_tmp)))
    for i in range(20):
        for j in range(20):
            matrix_final[i][j] = matrix_final[i][j]/matrix_tmp[j]
    tpc_vector = average(matrix_final, 1.0)
    return tpc_vector[0]


# Mixture of row and column transformation based feature encoders
def edp(input_matrix):
    """
    EDP feature encoder
    :param input_matrix:
    :return:
    """
    STEP = 2
    PART = 0
    ID = 1
    edp_matrix = [[0.0] * 20] * 1
    edp_matrix = np.array(edp_matrix)
    seq_cn = float(np.shape(input_matrix)[0])
    output_matrix = pre_handle_columns(input_matrix, STEP, PART, ID)
    output_matrix = np.array(output_matrix)
    for i in range(20):
        edp_matrix[0] = list(map(sum, zip(output_matrix[i], edp_matrix[0])))
    edp_matrix = np.divide(edp_matrix, (seq_cn-STEP)*20.0)
    return edp_matrix[0]


def rpssm(input_matrix):
    """
    RPSSM feature encoder
    :param input_matrix:
    :return:
    """
    rpssm_matrix = handle_mixed3(input_matrix)
    return rpssm_matrix[0]


def pse_pssm(input_matrix, alpha=1):
    """
    PSE PSSM feature encoder
    :param input_matrix:
    :param alpha:
    :return:
    """
    pse_pssm_matrix = handle_mixed(input_matrix, alpha)
    return pse_pssm_matrix[0]


def dp_pssm(input_matrix, alpha=5):
    """
    DP PSSM feature encoder
    :param input_matrix:
    :param alpha:
    :return:
    """
    dp_pssm_matrix = handle_mixed2(input_matrix, alpha)
    return dp_pssm_matrix[0]


def pssm_ac(input_matrix, group=10):
    """
    PSSM AC feature encoder
    :param input_matrix:
    :param group:
    :return:
    """
    ID = 0
    pssm_ac_matrix = correlation(input_matrix, ID, group)
    return pssm_ac_matrix[0]


def pssm_cc(input_matrix, group=10):
    """
    PSSM CC feature encoder
    :param input_matrix:
    :param group:
    :return:
    """
    ID = 1
    pssm_cc_matrix = correlation(input_matrix, ID, group)
    return pssm_cc_matrix[0]


# Feature encoding algorithms that combine the above encoders
def aadp_pssm(input_matrix):
    """
    aac and dpc combined encoder
    :param input_matrix:
    :return:
    """
    aac_pssm_matrix = aac_pssm(input_matrix)
    dpc_pssm_matrix = dpc_pssm(input_matrix)
    aac_pssm_matrix = np.array(aac_pssm_matrix)
    dpc_pssm_matrix = np.array(dpc_pssm_matrix)
    aadp_pssm_matrix = np.hstack((aac_pssm_matrix, dpc_pssm_matrix))
    return aadp_pssm_matrix


def aatp(input_matrix):
    """
    aac and tpc combined encoder
    :param input_matrix:
    :return:
    """
    aac_pssm_matrix = aac_pssm(input_matrix)
    tpc_matrix = tpc(input_matrix)
    aac_pssm_matrix = np.array(aac_pssm_matrix)
    tpc_matrix = np.array(tpc_matrix)
    aatp_matrix = np.hstack((aac_pssm_matrix, tpc_matrix))
    return aatp_matrix


def medp(input_matrix):
    """
    edp and eedp combined encoder
    :param input_matrix:
    :return:
    """
    edp_matrix = edp(input_matrix)
    eedp_matrix = eedp(input_matrix)
    edp_matrix = np.array(edp_matrix)
    eedp_matrix = np.array(eedp_matrix)
    medp_matrix = np.hstack((edp_matrix, eedp_matrix))
    return medp_matrix


def get_feature(pssm_dir, algo_type="aac_pssm", store_dir="./"):
    pssm_files = [os.path.join(pssm_dir, pf.name) for pf in os.scandir(pssm_dir) if pf.name.endswith(".pssm")]
    pssm_mat = list(map(read_pssm_matrix, pssm_files))
    features = np.array(list(map(eval(algo_type), pssm_mat)))

    protein_names = [re.sub(".pssm", "", pf.name) for pf in os.scandir(pssm_dir) if pf.name.endswith(".pssm")]

    if store_dir:
        output_dir_name = os.path.join(store_dir, "tmp/")
        os.makedirs(output_dir_name, exist_ok=True)
        with open(os.path.join(output_dir_name, algo_type+".csv"), "w") as f:
            for pn, feats in zip(protein_names, features):
                f.write(pn + ",")
                f.write(",".join(list(map(str, feats))))
                f.write("\n")
        return algo_type
    else:
        return zip(protein_names, features)


def get_all_features(pssm_dir, store_dir="./"):
    if store_dir:
        algo_types = [get_feature(pssm_dir, enc, store_dir) for enc in all_encoders]
        return algo_types
    else:
        return [get_feature(pssm_dir, enc, store_dir) for enc in all_encoders]


# create pssm profile function
def create_pssm_profile(seq_file, out_dir, psiblast_exec, database_prefix, num_threads=24):
    """
    A function to create psiblast or pssm profile for protein sequences
    :param seq_file: A csv file with name of the protein followed by its sequence separated by a comma
    :param out_dir: The directory where the user would like to store the pssm profiles of all the sequences
    :param psiblast_exec: The path of the psiblast executable. psiblast program needs to be installed
    :param database_prefix: The path of the indexed blast database directory prefix
    :param num_threads: Number of threads to use while creating the psiblast profile
    :return: The output directory where the psiblast/pssm profiles are stored
    """

    os.makedirs(out_dir, exist_ok=True)

    for line in open(seq_file).readlines():
        parsed_line = line.strip().split(",")
        name, sequence = parsed_line[0], parsed_line[1]
        with open(name + '.txt', 'w') as f:
            f.write('>' + name + '\n' + sequence + '\n')
        myCmd = f"{psiblast_exec} -query {name}.txt -db {database_prefix} -num_iterations 3 -num_threads " \
                f"{num_threads} -out_ascii_pssm {os.path.join(out_dir, name + '.pssm')}"
        print("Generating psiblast profile for protein: " + name)
        os.system(myCmd)
        os.remove(name + '.txt')
    return out_dir


if __name__ == "__main__":
    zip_list = get_all_features("./data/pssm_files", store_dir="")
    print(list(zip_list[0]))
