import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import sklearn.mixture as mixture
import math
import networkx as nx

# import kmapper as km
from sklearn.cluster import DBSCAN

import gudhi as gd
import statmapper as stm
from sklearn_tda import MapperComplex

def fpkm2tpm(fpkm, skip=0):
    """ Converts FPKM values to TPM values. TPM values allow better
    comparison of gene expression between individuals.

    Parameters
    ----------
    fpkm : data frame
        Data frame containing gene expression values.
        Rows must be genes while columns are different individuals

    skip : scalar
        Number of metadata columns to skip from the beginning of the dataframe

    Returns
    -------
    tpm : data frame
        A copy of the original dataframe with TPM values in place
    """

    total_exp = fpkm.iloc[:, skip:].sum(axis='index')
    tpm = fpkm.iloc[:, skip:]*1e6
    tpm = tpm.div(other=total_exp, axis='columns')

    return tpm

def _z_score_sample(sample, cutoff=2**-15, kernel='scott', num=1024):
    """
    Do NOT use it directly but with the z_score function.

    Parameters
    ----------
    sample : array-like
        Expression levels of all genes of a fixed cultivar/individual

    cutoff : scalar
        Genes whose expression is lower than the cutoff value will be
        considered as non-expressed

    kernel : str, callable
        bandwidth method for the Gaussian KDE. Look scipy.stats.gaussian_kde
        for more documentation.

    num : int
        Resolution to distinguish the mean of the KDE distribution

    Returns
    -------
    zFPKM : array-like
        Array with the z-scores of the sample

    mu : float
        Mean value from the fitted Gaussian distribution

    U : float
        Mean value of the expression levels above the 50% quantile

    sigma : float
        STD of the fitted Gaussian distribution

    """
    sample = sample[sample > cutoff]
    sample = np.log2(sample)
    kde = stats.gaussian_kde(sample)
    x = np.linspace(start=sample.min(), stop=sample.max(), num=num)
    y = kde(x)
    mu = x[np.argmax(y)]

    U = np.mean(sample[sample > mu])

    sigma = (U - mu) * math.sqrt(0.5*math.pi)

    zFPKM = (sample - mu) / sigma

    return zFPKM, mu, U, sigma

def z_score(df, skip=0, cutoff=-15, kernel='scott', num=1024):
    """Computes the z-scores of a given dataframe with TPM or FPKM values.
    This is simply to arrange all the TPM values in a histogram (fit a
    kernel-density estimate using Gaussian kernels to be more precise)
    and then fit a normal curve. The z-score simply how far the TPM value
    is from the mean. This is a good alternative if the plan is to zoom
    in the highest expressed genes for each individual. z-scores are unitless.

    Parameters
    ----------
    df : dataframe
        Data frame containing gene expression values.
        Rows must be genes while columns are different individuals

    skip : scalar
        Number of metadata columns to skip from the beginning of the dataframe

    cutoff : scalar
        Genes with expression less than 2^cutoff will be considered nonexpressed

    kernel : str, callable
        bandwidth method for the Gaussian KDE. Look scipy.stats.gaussian_kde
        for more documentation.

    num : int
        Resolution to distinguish the mean of the KDE distribution

    Returns
    -------
    zcopy : dataframe
        A copy of the original dataframe with z-scores in place
    """

    zcopy = df.copy()
    for col in zcopy.columns[skip:]:
        z,mu,U,sigma = _z_score_sample(df[col], cutoff=2**cutoff, kernel=kernel, num=num)
        zcopy[col] = cutoff
        zcopy.loc[z.index, col] = z
    return zcopy

def explore_correlation(fpkm, method='pearson'):
    """Observe pairwise correlation-related information

    Parameters
    ----------
    fpkm : dataframe
        Data frame containing gene expression values.
        Rows must be genes while columns are different individuals

    method : {‘pearson’, ‘kendall’, ‘spearman’}
        Method of correlation. Refer to pandas.DataFrame.corr

    Returns
    -------
    correlation : numpy array
        Square matrix of pairwise patient correlation

    vals : list of 1D-arrays
        Arrays the length of number of patients:
        vals[0] = min(correlation);
        vals[1] = mean(correlation);
        vals[2] = max(correlation);
        vals[3] = max - min;
    """

    correlations = np.array(fpkm.corr(method='pearson'))

    vals = []
    vals.append(np.amin(correlations, axis=1))

    np.fill_diagonal(correlations, 0)

    vals.append(np.mean(correlations, axis=1))
    vals.append(np.amax(correlations, axis=1))
    vals.append(vals[2] - vals[0])

    np.fill_diagonal(correlations, 1)

    return correlations, vals

def make_mapper(df, df_name, vals, idx, val_names, dst = './', labels=None,
                label_note=' ', eps=1e6, n_cubes=200, overlap=0.5):
    """ Compute a mapper graph for the given data, filter and clutering parameters

    Parameters
    ----------
    df : data frame
        Data frame containing gene expression related data.
        Rows must be genes while columns are different individuals

    df_name : string
        Name of the data frame (to be used in the title of the graph
        and the mapper HTML filename)

    vals : list of arrays
        From `explore_correlation`

    idx : {0,1,2,3} Filter index
        0 (min corr), 1 (mean corr), 2 (max corr) 3 (diff corr)

    val_names : list of strings
        List with description of different values evaluated

    labels : array-like
        Label and metadata information to color the mapper graph

    label_note : string
        Info on labels to put in the mapper title

    dst : string
        path to folder to save the HTML file

    eps : scalar
        DBSCAN clustering parameter

    n_cubes : scalar
        number of intervals

    overlap : 0.5
        percentage of overlap between intervals

    Returns
    -------
    graph : KeplerMapper object
    """

    if labels is None:
        labels = np.array(vals[idx])
        print('using default color labels')

    mapper = km.KeplerMapper(verbose=1)
    graph = mapper.map(vals[idx], np.array(df).T, clusterer=DBSCAN(eps=eps),
                       cover=km.Cover(n_cubes=n_cubes, perc_overlap=overlap))

    filename = df_name
    filtername = val_names[idx]

    filename += '_' + filtername + '.html'

    mapper_title = ' '.join(df_name.split('_')).title()
    mapper_title += ' : ' + label_note + ' : Filter='
    mapper_title += ' '.join(filtername.split('_')).title()

    mapper.visualize(graph, path_html=dst+filename, title=mapper_title, color_function=labels)

    return graph

def tscores(x, gmm):
    means = gmm['means']
    sds = gmm['stdeviations']
    s_assng = gmm['weights']

    mu = np.sum(s_assng * means, axis = 1)
    sigma = np.zeros(len(mu))
    for i in range(s_assng.shape[1]):
        sigma = sigma + s_assng[:,i] *(sds[i]**2 + np.power(means[i] - mu,2))

    T1 = np.sum(s_assng/sds, axis = 1)*(x - mu)
    T2 = (1.0/np.sqrt(np.sum(s_assng*np.power(sds,2), axis = 1)))*(x - mu)
    T3 = 1.0/np.sqrt(sigma) * (x - mu)

    h_assng = (s_assng == np.max(s_assng, axis=1)[: , np.newaxis]).astype(float)

    mu = np.sum(h_assng * means, axis = 1)
    T0 = np.sum(h_assng/sds, axis = 1)*(x - mu)

    return T0,T1,T2,T3

def gmm_fit(data, x, n_components=2, cov_type='spherical'):

    gmm = dict()
    estimator = mixture.GaussianMixture(n_components=n_components, covariance_type=cov_type, random_state=42)
    estimator.fit(data.values.reshape(-1,1))

    responsibilities = estimator.predict_proba(x.reshape(-1, 1))
    pdf = np.exp(estimator.score_samples(x.reshape(-1, 1)))
    pdf_individual = responsibilities * pdf[:, np.newaxis]

    weights = pdf_individual/pdf[:, np.newaxis]
    means = estimator.means_.squeeze()
    stdeviations = np.sqrt(estimator.covariances_)

    gmm['pdf'] = pdf
    gmm['pdf_individual'] = pdf_individual
    gmm['weights'] = weights
    gmm['means'] = means
    gmm['stdeviations'] = stdeviations

    return gmm

def plot_tsummary(dst, sample, x, T, gmm, data_label = 'data_label', dpi=100, save_fig=False):

    means_ = gmm['means']
    covariances_ = gmm['stdeviations']
    pdf = gmm['pdf']
    pdf_individual = gmm['pdf_individual']

    cmean = ['mediumblue', 'fuchsia']
    ccovs = ['cornflowerblue', 'palevioletred']

    gauss = stats.norm(0,1)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    lfont = 18
    axes[0,0].hist(sample.values, 50, density=True, color='slategrey',
                   histtype='stepfilled', alpha=0.4, label=data_label)

    for i in range(len(means_)):
        axes[0,0].axvline(means_[i], c=cmean[i],
                          label = 'mu = {:.3}'.format(means_[i]), ls='-.')
        axes[0,0].axvline(means_[i]+covariances_[i], c=ccovs[i],
                          label='sd = {:.3}'.format(covariances_[i]), ls=':')
        axes[0,0].axvline(means_[i]-covariances_[i], c=ccovs[i], ls=':')

    axes[0,0].plot(x, pdf, '-k', lw=3, label='GMM')
    axes[0,0].plot(x, pdf_individual, '--k', lw=2, label='individual')

    axes[0,0].set_xlabel('Log2({})'.format(data_label), fontsize=lfont)
    axes[0,0].set_ylabel('Frequency', fontsize=lfont)

    for i in range(len(T)):
        axes[0,1].plot(x, T[i], label='T{}'.format(i))


    for i in range(len(means_)):
        axes[0,1].axvline(means_[i], c=cmean[i], label = 'mu', ls='-.')
        axes[0,1].axvline(means_[i]+covariances_[i], c=ccovs[i], label='sd', ls=':')
        axes[0,1].axvline(means_[i]-covariances_[i], c=ccovs[i], ls=':')

    axes[0,1].set_xlabel('Log2({})'.format(data_label), fontsize=lfont)
    axes[0,1].set_ylabel('T-score', fontsize=lfont)
    axes[0,1].grid()

    i,j=2,0
    foo = np.linspace(T[j].min(),T[j].max(),1024)
    axes[0,i].hist(T[j], 50, density=True, color='slategrey',
                   histtype='stepfilled', alpha=0.4, label='T{}'.format(j))
    axes[0,i].plot(foo,gauss.pdf(foo), c='indigo', label='N(0,1)', lw=3)
    axes[0,i].set_xlabel('T{} score'.format(j), fontsize=lfont)


    for i,j in enumerate([1,2,3]):
        foo = np.linspace(T[j].min(),T[j].max(),1024)
        axes[1,i].hist(T[j], 50, density=True, color='slategrey',
                       histtype='stepfilled', alpha=0.4, label='T{}'.format(j))
        axes[1,i].plot(foo,gauss.pdf(foo), c='indigo', label='N(0,1)', lw=3)
        axes[1,i].set_xlabel('T{} score'.format(j), fontsize=lfont)

    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            axes[i,j].legend(fontsize=14)
            axes[i,j].tick_params(labelsize=14)

    fig.suptitle(sample.name, fontsize=22)
    plt.tight_layout()

    if save_fig:
        plt.savefig(dst + 'tsummary_' + data_label +'-_-'+sample.name+'.png',
                    facecolor='white', transparent=False,
                    dpi=dpi, format='png', bbox_inches='tight')
        plt.close()

def gmm_separation(gmm):
    means = np.argsort(gmm['means'])
    splus = gmm['means'] + gmm['stdeviations']
    sminus = gmm['means'] - gmm['stdeviations']
    return sminus[means[1]] - splus[means[0]]

def gmm_tscores(sample, cutoff=2**-15,
                n_components=2, cov_type='spherical',
                data_label='data', print_plot=False,
                save_fig=False, dpi=100, dst='./'):

    #idx = sample.sort_values().index
    idx = pd.Int64Index(sample.argsort().values)
    sample = sample[sample > cutoff]
    sample = np.log2(sample)
    x = sample.sort_values().values

    gmm = gmm_fit(sample, x, n_components, cov_type)
    T = tscores(x, gmm)

    foo = np.log2(cutoff)*np.ones(len(idx) - len(x), dtype=sample.dtype)
    tscore = [np.zeros(len(idx)) for i in range(len(T))]

    for i in range(len(T)):
        bar = np.hstack((foo,T[i]))
        tscore[i][idx] = bar
        #tdf[i][sample.name] = tscore[i]

    if print_plot:
        plot_tsummary(dst, sample, x, T, gmm, data_label=data_label, dpi=dpi, save_fig=save_fig)

    return tscore, gmm

def distinct_subjects(mapper_info, iternodes):
    subjects = set()
    for node in iternodes:
        sbjcts = set(mapper_info[node]['indices'])
        subjects |= sbjcts

    return subjects

def score_component_positions(mapper_info, path, num_subjects):
    main_comp = np.full(num_subjects, fill_value=0, dtype=float)
    leaves = [path[0], path[-1]]
    a , b = mapper_info[leaves[0]]['patch'][0] , mapper_info[leaves[1]]['patch'][0]

    if a > b:
        path = np.flip(np.array(path))

    delta = 2./(len(path) - 1)

    subjects = set()
    for i in range(len(path)):
        sbjcts = set(mapper_info[path[i]]['indices'])
        diff_sbjcts = sbjcts - subjects
        main_comp[np.array(list(diff_sbjcts)).astype(int)] += (-1 + delta*i)
        subjects |= sbjcts

    excluded = set(range(len(main_comp))) - subjects
    main_comp[np.array(list(excluded))] = np.nan

    return main_comp

def signif_strand(mapper_info, comp, leaves):
    ends = itertools.combinations(leaves, 2)
    maxweight = 0

    for left, right in ends:
        path = nx.shortest_path(comp, source=left, target=right)
        weight = 0
        for node in path:
            weight += mapper_info[node]['size']

        if weight > maxweight:
            maxpath = path
            maxweight = weight

    total_sbj = distinct_subjects(mapper_info, comp.nodes())
    maxpath_sbj = distinct_subjects(mapper_info, maxpath)

    return maxpath, total_sbj - maxpath_sbj

def draw_mapper(G, mapper_info, params, filename='mapper_graph.png', savefig=False, verbose=False):
    sett = distinct_subjects(mapper_info, G.nodes)
    intersections = []

    for edge in G.edges():
        subj0 = set(mapper_info[edge[0]]['indices'])
        subj1 = set(mapper_info[edge[1]]['indices'])
        intersection = subj0 & subj1
        intersections.append(len(intersection))

    nx.draw(G, pos=nx.kamada_kawai_layout(G),
            node_color=[mapper_info[node]["colors"][0] for node in G.nodes()],
            node_size=[3*mapper_info[node]['size'] for node in G.nodes()],
            linewidths=1,
            width = 0.2*np.array(intersections),
            node_shape='o',
            font_size=12)

    if verbose:
        plt.title('resol = {} ; gains = {}\n{}\n\nNumber of nodes: {}\nTotal distinct subjects: {}'
                  .format(params['resolutions'][0], params['gains'][0], params['clustering'],len(G.nodes), len(sett)))
    else:
        plt.title('Total distinct subjects: {}'.format(len(sett)), fontsize=20)

    if savefig:
        plt.savefig(filename,
                    dpi=100, format='png', bbox_inches='tight')
        #plt.close()

def bootstrap_scores(data, params, N=100):
    filter_func = params['filters']
    num_pts = len(data)
    chosen = np.zeros(len(data), dtype=int)
    bootscores = np.zeros((len(data), N), dtype=np.float)
    pruned = np.zeros(len(data), dtype=int)

    for bootstrap_id in range(N):
        single_node = True
        while single_node:
            # Randomly select points
            idxs = np.sort(np.unique(np.random.choice(num_pts, size=num_pts, replace=True)))
            selected = np.full(len(data), fill_value=False, dtype=bool)
            selected[idxs] = True
            chosen[np.unique(idxs)] += 1
            Xboot = data[idxs,:]
            f_boot = [filter_func[i] for i in idxs]
            params_boot = {k: params[k] for k in params.keys()}
            params_boot["filters"] = params["filters"][idxs,:]
            params_boot["colors"] = params["colors"][idxs,:]

            Mboot = MapperComplex(**params_boot).fit(Xboot)
            Gboot = stm.mapper2networkx(Mboot)
            conn_comps = [Gboot.subgraph(c) for c in sorted(nx.connected_components(Gboot), key=len, reverse=True)]
            leaves = [x for x in conn_comps[0].nodes() if conn_comps[0].degree(x) == 1]
            if len(leaves) > 1:
                single_node = False

        #print(bootstrap_id, '\tlen(leaves)', len(leaves))
        maxpath, pruned_sbj = signif_strand(Mboot.node_info_, conn_comps[0], leaves)
        if len(pruned_sbj) > 0:
            pruned[np.array(list(pruned_sbj))] += 1

        tmp_comp = score_component_positions(Mboot.node_info_, maxpath, len(idxs))

        bootscores[ selected, bootstrap_id] = tmp_comp
        bootscores[~selected, bootstrap_id] = np.nan

    return bootscores, chosen, pruned

