######!/usr/bin/env python3
"""
Contains required modules to process simulations

@author: ulas isildak
@e-mail: isildak.ulas [at] gmail.com
"""

import os
import sys
import csv
import allel
import random
import numpy as np
import pandas as pd
from PIL import Image

from sklearn.decomposition import PCA
#from sklearn.preprocessing import Imputer
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import train_test_split


def read_msms(filename, NCHROMS, N):
    """
    Reads msms file to an haplotype matrix
    Parameters:
        filename: full path and name of the .txt MSMS file
        NCHROMS: number of samples(haploid individuals, or chromosoms
        N: length of the simulated sequence(bp)
    Output:
        Returns an haplotype array, and an array containing positions
    """
    file = open(filename).readlines()
    if len(file) == 0:
        raise Exception('The file {} is empty'.format(filename.split('/')[-1]))
    # look for the // character in the file
    pointer = file.index('//\n') + 3
    # Get positions
    pos = file[pointer - 1].split()
    del pos[0]
    pos = np.array(pos, dtype='float')
    pos = pos * N
    # Get the number of genomic positions(determined be the number or pointers)
    n_columns = len(list(file[pointer])) - 1
    # Intialize the empty croms matrix: of type: object
    #croms = np.empty((NCHROMS, n_columns), dtype=np.object)
    croms = np.empty((NCHROMS, n_columns), dtype=object)
    # Fill the matrix with the simulated data
    for j in range(NCHROMS):
        f = list(file[pointer + j])
        del f[-1]
        F = np.array(f)
        croms[j, :] = F
    croms = croms.astype(int)
    return croms, pos


def rearrange_neutral(croms, pos, length):
    """
    rearranges neutral simulations such that each simulation results in
    <length>bp in length with selected mutation at the center
    Parameters:
        croms: input haplotype matrix
        pos: array containing position information for croms
        length: desired length of the output(bp)
    Returns:
        a haplotype matrix that is <length_out> bp in length and target snp at center,
        and an array containing new positions
    """
    freqs = np.true_divide(np.sum(croms, axis=0), croms.shape[0])
    # positions of mutations within [0.4,0.6]
    poss = pos[np.logical_and(freqs > 0.4, freqs < 0.6)]
    # position of target mutation
    pos_mut = poss[len(poss) // 2]
    # upper and lower boundaries of the region that will be selected
    up_bound = pos_mut + length / 2
    low_bound = pos_mut - length / 2
    target_range = np.logical_and(pos > low_bound, pos < up_bound)
    pos_new = pos[target_range] + length / 2 - pos_mut
    croms_new = croms[:, target_range]
    return croms_new, pos_new


def sort_freq(im_matrix):
    """
    This function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted
    by genetic similarity.
    Parameters:
        im_matrix: Array containing sequence data
    Returns:
        Sorted array containing sequence data
    """
    # u: Sorted Unique arrays
    # index: Index of 'im_matrix' that corresponds to each unique array
    # count: The number of instances each unique array appears in the 'im_matrix'
    u, index, count = np.unique(im_matrix, return_index=True, return_counts=True, axis=0)
    # b: Intitialised matrix the size of the original im_matrix[where new sorted data will be stored]
    b = np.zeros((np.size(im_matrix, 0), np.size(im_matrix, 1)), dtype=int)
    # c: Frequency table of unique arrays and the number of times they appear in the original 'im_matrix'
    c = np.stack((index, count), axis=-1)
    # The next line sorts the frequency table based mergesort algorithm
    c = c[c[:, 1].argsort(kind='mergesort')]
    pointer = np.size(im_matrix, 0) - 1
    for j in range(np.size(c, 0)):
        for conta in range(c[j, 1]):
            b[pointer, :] = im_matrix[c[j, 0]]
            pointer -= 1
    return b


def sort_min_diff(im_matrix):
    """
    This function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted
    by genetic similarity. this problem is NP, so here we use a nearest neighbors approx.  it's not perfect,
    but it's fast and generally performs ok.
    Implemented from https://github.com/flag0010/pop_gen_cnn/blob/master/sort.min.diff.py#L1
    Parameters:
        im_matrix: haplotype matrix (np array)
    Returns:
        Sorted numpy array
    """
    mb = NearestNeighbors(len(im_matrix), metric='manhattan').fit(im_matrix)
    v = mb.kneighbors(im_matrix)
    smallest = np.argmin(v[0].sum(axis=1))
    return im_matrix[v[1][smallest]]


def order_data(im_matrix, pos, sort, method):
    """
    Sorts haplotype matrix
    Parameters:
        im_matrix: input haplotype matrix
        pos: position of target SNP. not required if method = t
        sort: sorting method. either
            gen_sim: based on genetic similarity, or
            freq: based on frequency
        method: either
            t: together. sort the whole array together
            s: seperate. sort two haplotype groups seperately
    Returns:
        sorted haplotype matrix
    """
    if method == "t":
        if sort == "gen_sim":
            croms = sort_min_diff(im_matrix)
        elif sort == "freq":
            croms = sort_freq(im_matrix)
        else:
            raise ValueError("sort must be either 'freq' or 'gen_sim'")
    elif method == "s":
        if not isinstance(pos, int):
            raise ValueError("Position of the target SNP must be an integer")
        index_1 = (im_matrix[:, pos] == 1).reshape(im_matrix.shape[0], )
        index_0 = (im_matrix[:, pos] == 0).reshape(im_matrix.shape[0], )
        croms_1 = im_matrix[index_1, :]
        croms_0 = im_matrix[index_0, :]
        if sort == "gen_sim":
            croms1 = sort_min_diff(croms_1)
            croms0 = sort_min_diff(croms_0)
        elif sort == "freq":
            croms1 = sort_freq(croms_1)
            croms0 = sort_freq(croms_0)
        else:
            raise ValueError("sort must be either 'freq' or 'gen_sim'")
        croms = np.concatenate((croms0, croms1), axis=0)
    return croms


def sim_to_matrix(filename, NCHROMS, N, N_NE, sort, method):
    """
    Generates ordered haplotype_matrix from simulation results(must be ms format in .txt)
    Parameters:
        filename: full path and name of the .txt MSMS file
        NCHROMS: number of samples(haploid individuals, or chromosoms)
        N: length of the simulated sequence(bp)
        sort: sorting method. either:
            gen_sim: based on genetic similarity
            freq: based on frequency
        method: sorting method. either:
            t: together. sort the whole array together
            s: seperate. sort two haplotype groups seperately
    Returns:
        ordered haplotype matrix, where columns are positions and rows are individuals
    """
    if filename.split("/")[-1].startswith("NE"):
        crom, pos = read_msms(filename, NCHROMS, N_NE)
        croms, positions = rearrange_neutral(crom, pos, N)
    else:
        croms, positions = read_msms(filename, NCHROMS, N)
    pos = np.where(np.abs(positions - N / 2) < 1)[0]
    if len(pos) == 0:
        print(filename)
        print(positions)
        raise IndexError("Target SNP not found")
    if len(pos) > 1:
        print("Target SNP found at multiple positions:")
        print(positions[pos])
        pos = np.array([pos[0]], dtype='int64')
        print(filename)
    pos = int(pos[0])
    sorted_croms = order_data(croms, pos, sort, method)
    return sorted_croms


def matrix_to_image(croms, n_row, n_col):
    """
    Generates image from sorted haplotype matrix and resizes image into (n_row, ncol)
    Parameters:
        croms: Haplotype array
        n_row: number of rows
        n_col: number of cols
    Returns:
        resized image
    """
    # Generate image
    all1 = np.ones(croms.shape)
    cromx = all1 - croms
    bw_croms_uint8 = np.uint8(cromx)
    bw_croms_im = Image.fromarray(bw_croms_uint8 * 255, mode='L')
    # Resize
    im_resized = bw_croms_im.resize((n_col, n_row))
    return im_resized


def sim_to_image(path_to_sim, path_to_image, SIM_FROM, SIM_TO, NCHROMS, N, N_NE, img_dim=(128, 128),
                 clss=("NE", "IS", "OD", "FD"), sort="freq", method="s"):
    """
    Converts MSMS simulation output files into images
    Parameters:
        path_to_sim: path to the folder containing simulation files
        path_to_image: path to a folder in which output images will be saved
        SIM_FROM: number of simulations (replicate), starting from
        SIM_TO: number of simulations, until
        NCHROMS: number of samples(haploid individuals, or chromosomes)
        N: length of the simulated sequence(bp) for selection scenarios
        N_NE: length of simulated sequence(bp) for neutral scenario
        img_dim: image dimensions (nrow, ncol)
        clss: a tuple of target classes-
            "NE": neutral
            "IS": incomplete sweep
            "OD": overdominance
            "FD": negative freq-dependent selection
        sort: sorting algorithm. either
            gen_sim: based on genetic similarity
            freq: based on frequency
        method: either
            t: together. sort the whole array together
            s: seperate. sort two haplotype groups seperately
    """
    if sys.version_info >= (3, 6):
        with os.scandir(path_to_sim) as fdir:
            for file in fdir:
                if file.name.startswith(tuple(clss)) and file.is_file() and int(file.name.replace(".txt", "").split("_")[-1]) in range(SIM_FROM, SIM_TO+1):
                    croms = sim_to_matrix(file.path, NCHROMS, N, N_NE, sort=sort, method=method)
                    im_resized = matrix_to_image(croms, n_row=img_dim[0], n_col=img_dim[1])
                    im_resized.save("{}{}.bmp".format(path_to_image, file.name.replace(".txt", "")))
    else:
        files = [file for file in os.scandir(path_to_sim)
                 if file.is_file()
                 if file.name.startswith(tuple(clss))
                 if int(file.name.replace(".txt", "").split("_")[-1]) in range(SIM_FROM, SIM_TO + 1)]
        for file in files:
            croms = sim_to_matrix(file.path, NCHROMS, N, N_NE, sort=sort, method=method)
            im_resized = matrix_to_image(croms, n_row=img_dim[0], n_col=img_dim[1])
            im_resized.save("{}{}.bmp".format(path_to_image, file.name.replace(".txt", "")))
    return 0


def calc_median_r2(g):
    """Calculates median LD r^2"""
    gn = g.to_n_alt(fill=-1)
    LDr = allel.rogers_huff_r(gn)
    LDr2 = LDr ** 2
    median_r2 = np.nanmedian(LDr2)
    return median_r2


def calc_kelly_zns(g, n_pos):
    """Calculates Kelly's Zns statistic"""
    gn = g.to_n_alt(fill=-1)
    LDr = allel.rogers_huff_r(gn)
    LDr2 = LDr ** 2
    kellyzn = (np.nansum(LDr2) * 2.0) / (n_pos * (n_pos - 1.0))
    return kellyzn


def calc_pi(croms):
    """Calculates pi"""
    dis1 = []
    for i in range(croms.shape[0]):
        d1 = []
        for j in range(i + 1, croms.shape[0]):
            d1.append(sum(croms[i, :] != croms[j, :]))
        dis1.append(sum(d1))
    pi_est1 = (sum(dis1) / ((croms.shape[0] * (croms.shape[0] - 1.0)) / 2.0))
    return pi_est1


def calc_faywu_h(croms):
    """Calculates Fay and Wu's H statistic"""
    n_sam1 = croms.shape[0]
    counts1 = croms.sum(axis=0)
    S_i1 = []
    for i in range(1, n_sam1):
        S_i1.append(sum(counts1 == i))
    i = range(1, n_sam1)
    n_i = np.subtract(n_sam1, i)
    thetaP1 = sum((n_i * i * S_i1 * 2) / (n_sam1 * (n_sam1 - 1.0)))
    thetaH1 = sum((2 * np.multiply(S_i1, np.power(i, 2))) / (n_sam1 * (n_sam1 - 1.0)))
    Hstat1 = thetaP1 - thetaH1
    return Hstat1


def calc_fuli_f_star(croms):
    """Calculates Fu and Li's D* statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, range(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(range(1, n_sam1), 2)))
    an1 = an + np.true_divide(1, n_sam1)
    vfs = (((2 * (n_sam1 ** 3.0) + 110.0 * (n_sam1 ** 2.0) - 255.0 * n_sam1 + 153) / (
            9 * (n_sam1 ** 2.0) * (n_sam1 - 1.0))) + ((2 * (n_sam1 - 1.0) * an) / (n_sam1 ** 2.0)) - (
                   (8.0 * bn) / n_sam1)) / ((an ** 2.0) + bn)
    ufs = ((n_sam1 / (n_sam1 + 1.0) + (n_sam1 + 1.0) / (3 * (n_sam1 - 1.0)) - 4.0 / (
            n_sam1 * (n_sam1 - 1.0)) + ((2 * (n_sam1 + 1.0)) / ((n_sam1 - 1.0) ** 2)) * (
                    an1 - ((2.0 * n_sam1) / (n_sam1 + 1.0)))) / an) - vfs
    pi_est = calc_pi(croms)
    ss = sum(np.sum(croms, axis=0) == 1)
    Fstar1 = (pi_est - (((n_sam1 - 1.0) / n_sam1) * ss)) / ((ufs * n_pos1 + vfs * (n_pos1 ** 2.0)) ** 0.5)
    return Fstar1


def calc_fuli_d_star(croms):
    """Calculates Fu and Li's D* statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, range(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(range(1, n_sam1), 2)))
    an1 = an + np.true_divide(1, n_sam1)
    cn = (2 * (((n_sam1 * an) - 2 * (n_sam1 - 1))) / ((n_sam1 - 1) * (n_sam1 - 2)))
    dn = (cn + np.true_divide((n_sam1 - 2), ((n_sam1 - 1) ** 2)) + np.true_divide(2, (n_sam1 - 1)) * (
            3.0 / 2 - (2 * an1 - 3) / (n_sam1 - 2) - 1.0 / n_sam1))
    vds = (((n_sam1 / (n_sam1 - 1.0)) ** 2) * bn + (an ** 2) * dn - 2 * (n_sam1 * an * (an + 1)) / (
            (n_sam1 - 1.0) ** 2)) / (an ** 2 + bn)
    uds = ((n_sam1 / (n_sam1 - 1.0)) * (an - n_sam1 / (n_sam1 - 1.0))) - vds
    ss = sum(np.sum(croms, axis=0) == 1)
    Dstar1 = ((n_sam1 / (n_sam1 - 1.0)) * n_pos1 - (an * ss)) / (uds * n_pos1 + vds * (n_pos1 ^ 2)) ** 0.5
    return Dstar1


def calc_zeng_e(croms):
    """Calculates Zeng et al's E statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, range(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(range(1, n_sam1), 2)))
    counts1 = croms.sum(axis=0)
    S_i1 = []
    for i in range(1, n_sam1):
        S_i1.append(sum(counts1 == i))
    thetaW = n_pos1 / an
    thetaL = np.sum(np.multiply(S_i1, range(1, n_sam1))) / (n_sam1 - 1.0)
    theta2 = (n_pos1 * (n_pos1 - 1.0)) / (an ** 2 + bn)
    var1 = (n_sam1 / (2.0 * (n_sam1 - 1.0)) - 1.0 / an) * thetaW
    var2 = theta2 * (bn / (an ** 2.0)) + 2 * bn * (n_sam1 / (n_sam1 - 1.0)) ** 2.0 - (
            2.0 * (n_sam1 * bn - n_sam1 + 1.0)) / ((n_sam1 - 1.0) * an) - (3.0 * n_sam1 + 1.0) / (
               (n_sam1 - 1.0))
    varlw = var1 + var2
    ZengE1 = (thetaL - thetaW) / (varlw) ** 0.5
    return ZengE1


def calc_rageddness(croms):
    """Calculates rageddness index"""
    mist = []
    for i in range(croms.shape[0] - 1):
        for j in range(i + 1, croms.shape[0]):
            mist.append(sum(croms[i, :] != croms[j, :]))
    mist = np.array(mist)
    lnt = mist.shape[0]
    fclass = []
    for i in range(1, np.max(mist) + 2):
        fclass.append((np.true_divide(sum(mist == i), lnt) - np.true_divide(sum(mist == (i - 1)), lnt)) ** 2)
    rgd1 = np.sum(fclass)
    return rgd1




def sim_to_stats(path_to_sim, path_to_stat, clss, NCHROMS, SIM_FROM, SIM_TO, N, N_NE):
    """
    Calculates summary statistics for simulation outputs. Creates .csv file at path_to_stat containing summary
    statistics for specified simulations.
    Parameters:
        path_to_sim: Path to directory where the simulation files exist
        path_to_stat: Path to directory where the summary statistics will be stored
        clss: Class of the simulation(either FD, OD, IS or NE):
            -FD: negative-frequency dependent selection
            -OD: over dominance
            -IS: incomplete sweep
            -NE: neutral
        NCHROMS: number of samples(haploid individuals, or chromosoms)
        SIM_FROM: number of simulations (replicate) -starting from
        SIM_TO: number of simulations -until
        N: length of the simulated sequence(bp) for selection scenarios
        N_NE: length of simulated sequence(bp) for neutral scenario
    """
    if SIM_FROM == 1:
        once = 0
    else:
        once = 1
    files = [file for file in os.scandir(path_to_sim)
             if file.is_file()
             if file.name.startswith(tuple(clss))
             if int(file.name.replace(".txt", "").split("_")[-1]) in range(SIM_FROM, SIM_TO + 1)]
    files = [f.path for f in files]
    for file in sorted(files):
        if file.split("/")[-1].startswith("NE"):
            crom, pos = read_msms(file, NCHROMS, N_NE)
            croms, positions = rearrange_neutral(crom, pos, N)
        else:
            croms, positions = read_msms(file, NCHROMS, N)
        fname = file.split("/")[-1].replace(".txt", "")
        labs, stats = sum_stats(croms, positions, NCHROMS, fname)
        f = open("{}{}.csv".format(path_to_stat, fname.split('_')[0]), 'a+')
        with f:
            writer = csv.DictWriter(f, fieldnames=labs)
            if once == 0:
                writer.writeheader()
                writer.writerow(dict(zip(labs, stats)))
                once = 1
            else:
                writer.writerow(dict(zip(labs, stats)))
    return 0






def sum_stats(croms, pos, NCHROMS, run, sampling_time):
    """
    Calculates summary statistics
    Parameters:
        croms: haplotype matrix
        pos: positions
        sname: simulation name
        NCHROMS: number of chromosomes
    Returns:
        A list of labels (names of statistics)
        A list of values
    """
    # SUMMARY STATISTICS
    # REGION 1: full region ([0bp:50000bp])
    pos1 = pos[np.logical_and(pos > 0, pos < 50000)]
    croms1 = croms[:, np.logical_and(pos > 0, pos < 50000).tolist()]
    n_pos1 = np.size(croms1, 1)
    if n_pos1 == 0:
        raise ValueError("Region 1 contains no positions. Check the input data.")
        
    freq1 = np.true_divide(croms1.sum(axis=0), NCHROMS)
    freq1 = np.array(freq1)
    haplos = np.transpose(croms1)
    h1 = allel.HaplotypeArray(haplos)
    ac1 = h1.count_alleles()
    g1 = h1.to_genotypes(ploidy=2, copy=True)
    # mean_pairwise_distance
    mean_mean_pwise_dis1 = np.mean(allel.mean_pairwise_difference(ac1))
    median_mean_pwise_dis1 = np.median(allel.mean_pairwise_difference(ac1))
    max_mean_pwise_dis1 = np.max(allel.mean_pairwise_difference(ac1))
    # tajimasd
    TjD1 = allel.tajima_d(ac1)
    # watterson
    theta_hat_w1 = allel.watterson_theta(pos1, ac1)
    # heterogeneity
    obs_het1 = allel.heterozygosity_observed(g1)
    af1 = ac1.to_frequencies()
    exp_het1 = allel.heterozygosity_expected(af1, ploidy=2)
    mean_obs_het1 = np.mean(obs_het1)
    median_obs_het1 = np.median(obs_het1)
    max_obs_het1 = np.max(obs_het1)
    # return 0, if 0/0 encountered
    ob_exp_het1 = np.divide(obs_het1, exp_het1, out=np.zeros_like(obs_het1), where=exp_het1 != 0)
    mean_obs_exp1 = np.nanmean(ob_exp_het1)
    median_obs_exp1 = np.nanmedian(ob_exp_het1)
    max_obs_exp1 = np.nanmax(ob_exp_het1)
    # LD r
    median_r21 = calc_median_r2(g1)
    # Haplotype_stats
    hh1 = allel.garud_h(h1)
    h11 = hh1[0]
    h121 = hh1[1]
    h1231 = hh1[2]
    h2_h11 = hh1[3]
    n_hap1 = np.unique(croms1, axis=0).shape[0]
    hap_div1 = allel.haplotype_diversity(h1)
    ehh1 = allel.ehh_decay(h1)
    mean_ehh1 = np.mean(ehh1)
    median_ehh1 = np.median(ehh1)
    ihs1 = allel.ihs(h1, pos1, include_edges=True)
    median_ihs1 = np.nanmedian(ihs1)
    # nsl
    nsl1 = allel.nsl(h1)
    max_nsl1 = np.nanmax(nsl1)
    median_nsl1 = np.nanmedian(nsl1)
    # NCD
    tf = 0.5
    freq11 = freq1[freq1 < 1]
    n1 = freq11.shape[0]
    ncd11 = (sum((freq11 - tf) ** 2) / n1) ** 0.5
    # kellyZns
    kellyzn1 = calc_kelly_zns(g1, n_pos1)
    # pi
    pi_est1 = calc_pi(croms1)
    # FayWusH
    Hstat1 = calc_faywu_h(croms1)
    # of singletons
    Ss1 = sum(np.sum(croms1, axis=0) == 1)
    # fu_li Dstar
    Dstar1 = calc_fuli_d_star(croms1)
    # fu_li Fstar
    Fstar1 = calc_fuli_f_star(croms1)
    # Zeng_E
    ZengE1 = calc_zeng_e(croms1)
    # rageddness index
    rgd1 = calc_rageddness(croms1)

    stats = ['Neutral', str(run), str(sampling_time), mean_mean_pwise_dis1, median_mean_pwise_dis1, max_mean_pwise_dis1,
             TjD1, theta_hat_w1, mean_obs_het1, median_obs_het1, max_obs_het1, mean_obs_exp1, median_obs_exp1,
             max_obs_exp1, median_r21, h11, h121, h1231, h2_h11, hap_div1, n_hap1, mean_ehh1, median_ehh1, median_ihs1,
             max_nsl1, median_nsl1, ncd11, kellyzn1, pi_est1, Hstat1, Ss1, Dstar1, Fstar1, ZengE1, rgd1]

    labs = ['Selection_mode', 'H_S', 'Run', 'Generation', 'Mean(MeanPwiseDist)1', 'Median(MeanPwiseDist)1', 'Max(MeanPwiseDist)1',
            'Tajimas D1', 'Watterson1', 'Mean(ObservedHet)1', 'Median(ObservedHet)1', 'Max(ObservedHet)1',
            'Mean(Obs/Exp Het)1', 'Median(Obs/Exp Het)1', 'Max(Obs/Exp Het)1', 'Median(r2)1', 'H1_1', 'H12_1',
            'H123_1', 'H2/H1_1', 'Haplotype Diversity1', '# of Hap1', 'Mean(EHH)1', 'Median(EHH)1', 'Median(ihs)1',
            'Max(nsl)1', 'Median(nsl)1', 'NCD1_1', 'KellyZns1', 'pi1', 'faywuH1', '#ofSingletons1',
            'Dstar1', 'Fstar1', 'ZengE1', 'Rageddnes1']

    return labs, stats



FILENAME = sys.argv[1]
RUN = sys.argv[2]
SAMPLING_TIME = sys.argv[3]

print(FILENAME, RUN, SAMPLING_TIME)

# - 16 haploid chromosome (8 diploid individuals since 8*2=16)
# - 50000 base pair genomic region length
chroms, positions = read_msms(FILENAME, 16, 50000)


#statistics = sum_stats(chroms, positions, 2000, 'NA_NA_' + RUN)[1]
statistics = sum_stats(chroms, positions, 16, RUN, SAMPLING_TIME)[1]

outf = open('/mnt/NEOGENE1/projects/selection_2023/imputation/imputed_dataset/snp_by_snp2/P1/filtered_vcf_files/sum_stat_deneme.txt', 'a')
print('\t'.join(str(element) for element in statistics), file = outf)
outf.close()


####################################################################
####################################################################

