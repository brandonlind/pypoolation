"""
A python version of some scripts from popoolation.

TODO: if Tajima's D is selected, return pi and theta (and maybe ddivisor)
"""

import os, sys, argparse, shutil, subprocess, pandas as pd, threading, ipyparallel, time
import pickle
from os import path as op
from typing import Union


def read(file:str, lines=True) -> Union[str, list]:
    """Read contents of file."""
    with open(file, 'r') as o:
        text = o.read()

    if lines is True:
        ret = text.split("\n")
    else:
        ret = text

    return ret


def pklload(path:str):
    """Load object from a .pkl file."""
    pkl = pickle.load(open(path, 'rb'))
    return pkl


def pkldump(obj, f:str) -> None:
    """Save object to .pkl file."""
    with open(f, 'wb') as o:
        pickle.dump(obj, o, protocol=pickle.HIGHEST_PROTOCOL)


def get_client(profile='default') -> tuple:
    """Get lview,dview from ipcluster."""
    rc = ipyparallel.Client(profile=profile)
    dview = rc[:]
    lview = rc.load_balanced_view()

    return lview, dview


def watch_async(jobs:list, phase=None) -> None:
    """Wait until jobs are done executing, show progress bar."""
    from tqdm import trange

    print(ColorText(f"\nWatching {len(jobs)} {phase} jobs ...").bold())
    
    job_idx = list(range(len(jobs)))
    for i in trange(len(jobs)):
        count = 0
        while count < (i+1):
            count = len(jobs) - len(job_idx)
            for j in job_idx:
                if jobs[j].ready():
                    count += 1
                    job_idx.remove(j)


class ColorText:
    """
    Use ANSI escape sequences to print colors +/- bold/underline to bash terminal.
    """
    def __init__(self, text:str):
        self.text = text
        self.ending = '\033[0m'
        self.colors = []

    def __str__(self):
        return self.text

    def bold(self):
        self.text = '\033[1m' + self.text + self.ending
        return self

    def underline(self):
        self.text = '\033[4m' + self.text + self.ending
        return self

    def green(self):
        self.text = '\033[92m' + self.text + self.ending
        self.colors.append('green')
        return self

    def purple(self):
        self.text = '\033[95m' + self.text + self.ending
        self.colors.append('purple')
        return self

    def blue(self):
        self.text = '\033[94m' + self.text + self.ending
        self.colors.append('blue')
        return self

    def warn(self):
        self.text = '\033[93m' + self.text + self.ending
        self.colors.append('yellow')
        return self

    def fail(self):
        self.text = '\033[91m' + self.text + self.ending
        self.colors.append('red')
        return self


def get_pars():
    """
    Parse input flags.

    # TODO check arg descriptions, and if they're actually used.
    """
    parser = argparse.ArgumentParser(description=print(mytext),
                                     add_help=True,
                                     formatter_class=argparse.RawTextHelpFormatter)
    requiredNAMED = parser.add_argument_group('required arguments')
    requiredNAMED.add_argument("-i", "--input",
                               required=True,
                               default=None,
                               dest="input",
                               type=str,
                               help='''/path/to/VariantsToTable_output.txt
It is assumed that there is either a 'CHROM' or 'unstitched_chrom'
column (ref.fa record name), and either a 'locus' or 'unstitched_locus' column.
The 'locus' column elements are the hyphen-separated
CHROM-POS. If the 'unstitched_chrom' column is present, the code will use the
'unstitched_locus' column for SNP names, otherwise 'CHROM' and 'locus'. The
'unstitched_locus' elements are therefore the hyphen-separated 
unstitched_locus-unstitched_pos. DP, AD, and RD columns from VarScan are also 
assumed.
''')
    requiredNAMED.add_argument("-o","--outdir",
                               required=True,
                               default=None,
                               dest="outdir",
                               type=str,
                               help="/path/to/pypoolation_output_dir/")
    requiredNAMED.add_argument("-m","--measure",
                               required=True,
                               default=None,
                               dest="measure",
                               type=str,
                               help='''pi OR theta OR D.''')
    requiredNAMED.add_argument("-p","--ploidy",
                               required=True,
                               default=None,
                               dest="ploidyfile",
                               type=str,
                               help='''/path/to/the/ploidy.pkl file output by the VarScan pipeline. This is a python
dictionary with key=pool_name, value=dict with key=pop, value=ploidy. The code
will prompt for pool_name. pop(s) can be input with --whichpops.''')
    requiredNAMED.add_argument("-e","--engines",
                               required=True,
                               default=None,
                               dest="engines",
                               type=int,
                               help="The number of ipcluster engines that will be launched.")
    parser.add_argument("--ipcluster-profile",
                        required=False,
                        default='default',
                        dest="profile",
                        type=str,
                        help="The ipcluster profile name with which to start engines. Default: 'default'")
    parser.add_argument("--which-pops",
                        required=False,
                        nargs='+',
                        default=None,
                        dest="whichpops",
                        type=str,
                        help='''Specific populations you would like statistics calculated \
- entered as a space separated list. Default=all''')
    parser.add_argument("--min-count",
                        required=False,
                        default=2,
                        dest="mincount",
                        type=int,
                        help="The minimum count of the minor allele. Default=2")
    parser.add_argument("--min-coverage",
                        required=False,
                        default=8,
                        dest="mincov",
                        type=int,
                        help='''The minimum coverage of a site. Sites with a lower coverage
will not be considered. Default=8''')
    parser.add_argument("--max-coverage",
                        required=False,
                        default=1000,
                        dest="maxcov",
                        type=int,
                        help='''The maximum coverage of a site. Sites with greater values
will not be considered. Default=1000''')
    parser.add_argument("--min-coverage-fraction",
                        required=False,
                        default=0.0,
                        dest="mincovfraction",
                        type=float,
                        help='''The minimum fraction of a window having SNPs between min-coverage
and max-coverage (fraction = num snps / windowsize). For the
CoAdapTree data, this strategy was most computationally
efficient. Data was further filtered after output.  Default=0.0''')
    parser.add_argument("--window-size",
                        required=False,
                        default=1000,
                        dest="windowsize",
                        type=int,
                        help='''The size of the sliding window. Windows are centered on
individual SNPs, so any SNP within 0.5*windowsize bps
from the current SNP will be included in the window. Default=1000bp''')
    parser.add_argument("--chromfile",
                        required=False,
                        default=None,
                        dest="chromfile",
                        type=str,
                        help='''Path to .txt file with each line a chromosome/contig name
that is wished to be evaluated. All SNPs on this 
chromosome/contig will be evaluated. It is assumed that
these chromosomes are in the "CHROM" or "unstitched_chrom"
column of the input file. If 'unstitched_chrom' column
is present in input file, then this is the default column
to use to parallelize input, otherwise CHROM is used.
Default to analyze=All.''')

    # check flags
    args = parser.parse_args()
    if not op.exists(args.outdir):
        print(ColorText(f"FAIL: the directory for the output file(s) does not exist.").fail())
        print(ColorText(f"FAIL: please create this directory: %s" % args.outdir).fail())
        print(ColorText("exiting pypoolation.py").fail())
        exit()

    # make sure input and ploidyfile exist
    nopath = []
    for x in [args.input, args.ploidyfile]:  # TODO: check for $HOME or other bash vars in path
        if not op.exists(x):
            nopath.append(x)

    # if input or ploidy file do not exist:
    if len(nopath) > 0:
        print(ColorText("FAIL: The following path%s do not exist:" % "s" if len(nopath) > 1 else "").fail())
        [print(ColorText("\nFAIL: %s" % f).fail()) for f in nopath]
    
    print('args = ', args)

    return args


def askforinput(msg='Do you want to proceed?', tab='', newline='\n'):
    """Ask for input; if msg is default and input is no, exit."""
    while True:
        inp = input(ColorText(f"{newline}{tab}INPUT NEEDED: {msg} \n{tab}(yes | no): ").warn().__str__()).lower()
        if inp in ['yes', 'no']:
            if inp == 'no' and msg=='Do you want to proceed?':
                print(ColorText('exiting %s' % sys.argv[0]).fail())
                exit()
            break
        else:
            print(ColorText("Please respond with 'yes' or 'no'").fail())

    return inp


def get_datatable(args, pop, chroms=None) -> (pd.DataFrame, str):
    """Load --input datatable.txt."""
    print(ColorText(f"\nReading in SNP datatable for {pop}...").bold())

    # first get column names so that large dataframes are read in faster
    allcols = subprocess.check_output([shutil.which('head'), '-1', args.input]).decode('utf-8').split('\n')[0].split('\t')
    # get columns that are needed by pypoolation
    cols = [col for col in allcols if any(['chrom' in col.lower(),
                                           'pos' in col.lower(),
                                           f'{pop}.' in col,
                                           'locus' in col.lower()])]
    # use these cols to read in dataframe
    df = pd.read_table(args.input, usecols=cols)

    # see if user supplied specific chromosomes/contigs to analyze
    if args.chromfile is not None:
        if op.exists(args.chromfile):
            chroms = read(args.chromfile)
        else:
            print(ColorText(f"\nFAIL: chromsome file does not exist: {args.chromfile}").fail())
            print(ColorText("exiting pypoolation.py").fail())
            exit()

    # determine which column to use for parallelization
    chromcol = 'unstitched_chrom' if 'unstitched_chrom' in df.columns else 'CHROM'
    if chromcol not in df.columns:
        print(ColorText(f"FAIL: {chromcol} is not in the input columns").fail())
        print(ColorText("exiting pypoolation.py").fail.bold())
        exit()
    print(f'\tUsing the {chromcol} column to parallelize data')
    
    # reduce df to pop that matters
    keepcols = [col for col in df.columns if '.' not in col or pop in col]
    df = df[keepcols].copy()

    # reduce rows by removing missing data (np.nans)
    df = df[df[f'{pop}.FREQ']==df[f'{pop}.FREQ']].copy()  # only keep filtered SNPs (filtered SNPs don't have a FREQ)

    # set index to chromcol for db-like lookups
    df.index = df[chromcol].tolist()

    # reduce df to only chroms that user specified
    if chroms is not None:
        # make sure we're using the right column
        if len(set(df[chromcol]).intersection(chroms)) == 0:
            print(ColorText(f"\t\tWARN: None of the chromosomes/contigs in the chromosome file were found in the input file's {chromcol} column. Choosing to proceed will result in all chromosomes/contigs being analyzed.").warn())
            askforinput(tab='\t\t', newline='')
        else:
            df = df[df[chromcol].isin(chroms)].copy()

    return df, chromcol


def wait_for_engines(engines, profile):
    """Reload engines until number matches input engines arg."""
    lview = []
    dview = []
    count = 1
    while any([len(lview) != engines, len(dview) != engines]):
        if count % 30 == 0:
            # if waiting too long..
            # TODO: if found engines = 0, no reason to ask, if they continue it will fail
            print('count = ', count)
            print(ColorText("\tFAIL: Waited too long for engines.").fail())
            print(ColorText("\tFAIL: Make sure that if any cluster is running, the -e arg matches the number of engines.").fail())
            print(ColorText("\tFAIL: In some cases, not all expected engines can start on a busy server.").fail())
            print(ColorText("\tFAIL: Therefore, it may be the case that available engines will be less than requested.").fail())
            print(ColorText("\tFAIL: pypoolation found %s engines, with -e set to %s" % (len(lview), engines)).fail())
            answer = askforinput(msg='Would you like to continue with %s engines? (choosing no will wait another 60 seconds)' % len(lview), tab='\t', newline='')
            if answer == 'yes':
                break
        try:
            lview,dview = get_client(profile=profile)
        except (OSError, ipyparallel.error.NoEnginesRegistered, ipyparallel.error.TimeoutError):
            lview = []
            dview = []
        time.sleep(2)
        count += 1

    print('\tReturning lview,dview (%s engines) ...' % len(lview))

    return lview,dview


def launch_engines(engines, profile):
    """Launch ipcluster with engines under profile."""
    print(ColorText(f"\nLaunching ipcluster with {engines} engines...").bold())

    def _launch(engines, profile):
        subprocess.call([shutil.which('ipcluster'), 'start', '-n', str(engines), '--daemonize'])

    # first see if a cluster has already been started
    started = False
    try:
        print("\tLooking for existing engines ...")
        lview,dview = get_client(profile=profile)
        if len(lview) != engines:
            lview,dview = wait_for_engines(engines, profile)
        started = True
    except (OSError, ipyparallel.error.NoEnginesRegistered, ipyparallel.error.TimeoutError):
        print("\tNo engines found ...")

    # if not, launch 'em
    if started is False:
        print("\tLaunching engines ...")
        # pid = subprocess.Popen([shutil.which('ipcluster'), 'start', '-n', str(engines)]).pid
        x = threading.Thread(target=_launch, args=(engines,profile,), daemon=True)
        x.daemon=True
        x.start()
        lview,dview = wait_for_engines(engines, profile)

    return lview,dview


def attach_data(**kwargs) -> None:
    """Load object to engines."""
    import time

    num_engines = len(kwargs['dview'])
    print(ColorText("\nAdding data to engines ...").bold())
    print(ColorText("\tWARN: Watch available mem in another terminal window: 'watch free -h'").warn())
    print(ColorText("\tWARN: If available mem gets too low, kill engines and restart pypoolation with fewer engines: 'ipcluster stop'").warn())
    for key,value in kwargs.items():
        if key != 'dview':
            print(f'\tLoading {key} ({value.__class__.__name__}) to {num_engines} engines')
            kwargs['dview'][key] = value
            time.sleep(1)
    time.sleep(10)

    return None


def uni(lst) -> list:
    """Return unqiue items from lst."""
    return list(set(lst))


# def get_combos(pop, df, ploidy, dview):
#     """Get unique combinations of minor_count (b), pop ploidy (n), and coverage (M)."""
#     # get unique mincounts (b)  # TODO: don't worry about anthing less than args.mincount
#     df[f'{pop}.minor_count'] = np.nan
#     minors = df[f'{pop}.AD'] < df[f'{pop}.RD']
#     df.loc[minors, f'{pop}.minor_count'] = df.loc[minors, f'{pop}.AD']
#     df.loc[~minors, f'{pop}.minor_count'] = df.loc[~minors, f'{pop}.RD']
#     uni_b = uni([x for x in df[f'{pop}.minor_count'].astype(int) if x==x and x!=0])
#     print('uni_b = ', uni_b)

#     # get unique poolsizes (from ploidy)
#     uni_n = uni(ploidy.values())

#     # get unique coverages
#     uni_M = uni([x for x in df[f'{pop}.DP'].astype(int) if x==x])

#     # get unique combinations
#     uni_combos = list(itertools.product(range(len(uni_b)),
#                                         range(len(uni_M)),
#                                         range(len(uni_n))))

#     # attache data to engines
#     attach_data(uni_b=uni_b, uni_M=uni_M, uni_n=uni_n, dview=dview)

#     return uni_combos, uni_b, uni_M, uni_n


# def pidiv_iterator(chunk, uni_b, uni_M, uni_n):
#     """For each chunk (a unique combo of b, M, and n), calculate pidiv_buffer."""
#     import varmath

#     pi_buffer = {}
#     for b_idx, M_idx, n_idx in chunk:
#         b = int(uni_b[int(b_idx)])
#         M = int(uni_M[int(M_idx)])
#         n = int(uni_n[int(n_idx)])
#         key = f"{b}-{n}-{M}"
#         div, buffdict = varmath.get_pidiv_buffer(b, n, M)
#         pi_buffer[key] = div

#     return pi_buffer


# def send_get_pidiv_buffer(**kwargs):
#     """Parallelize varmath.get_pidiv_buffer()."""
#     import math

#     # get unique combos of b, n, and M
#     uni_combos, uni_b, uni_M, uni_n = get_combos(kwargs['pop'], kwargs['snps'], kwargs['ploidy'], kwargs['dview'])
#     for name,lst in zip(['uni_combos', 'uni_b', 'uni_M', 'uni_n'],[uni_combos, uni_b, uni_M, uni_n]):
#         print(name, len(lst))
#     print(uni_combos[:5])
    
#     # send out jobs
#     maxjobs = math.ceil(len(uni_combos)/kwargs['engines'])
#     print('maxjobs = ', maxjobs)
#     jobs = []
#     chunk = []
#     for i,combo in enumerate(uni_combos):
#         chunk.append(combo)
#         if len(chunk) == maxjobs or (i+1) == len(uni_combos):
#             jobs.append(kwargs['lview'].apply_async(pidiv_iterator, *(chunk, uni_b, uni_M, uni_n)))
# #             jobs.append(pidiv_iterator(chunk, uni_b, uni_M, uni_n))
#             chunk = []

#     # wait for jobs to finish
#     watch_async(jobs, 'pidiv_iterator')

#     # get returns
#     buffdict = {}
#     for j in jobs:
#         print('len(j.r) = ', len(j.r))
#         print('type(j.r) = ', type(j.r))
#         buffdict.update(j.r)
#     print('buffdict = ', list(buffdict.keys())[:5])

#     return buffdict


def choose_pool(ploidy:dict, args, pool=None) -> dict:
    """Choose which the pool to use as a key to the ploidy dict."""
    keys = list(ploidy.keys())
    if len(keys) == 1:
        # return the value of the dict using the only key
        return ploidy[keys[0]]

    print(ColorText('\nPlease choose a pool that contains the population of interest.').bold())
    nums = []
    for i,pool_name in enumerate(keys):
        print('\t%s %s' % (i, pool_name))
        nums.append(i)

    while True:
        inp = int(input(ColorText("\tINPUT NEEDED: Choose file by number: ").warn()).lower())
        if inp in nums:
            pool = keys[inp]
            break
        else:
            print(ColorText("\tPlease respond with a number from above.").fail())

    # make sure they've chosen at least one account
    while pool is None:
        print(ColorText("\tFAIL: You need to specify at least one pool. Revisiting options...").fail())
        pool = choose_pool(ploidy, args, keep=None)

    return ploidy[pool]


def get_ploidy(args) -> dict:
    """Get the ploidy of the populations of interest, reduce ploidy pkl."""
    # have user choose key to dict
    ploidy = choose_pool(pklload(args.ploidyfile), args)

    # keep only the pops from --which-pops flag
    for pop in list(ploidy.keys()):
        if not pop in args.whichpops:
            ploidy.pop(pop)

    # make sure outcome makes sense
    if not len(ploidy.keys()) == len(args.whichpops):
        print(ColorText("FAIL: The populations that were provided were not all found in ploidy.pkl").fail())
        print(ColorText("FAIL: Here are the populations that were provided that were not found:").fail())
        for pop in set(args.whichpops) - set(ploidy.keys()):
            print(ColorText(f"\t{pop}").fail())

    return ploidy


def get_windows(chrom, **kwargs) -> dict:
    """
    Get info for all SNPs meeting criteria for a given window size on a specific chrom.

    A window is centered on a SNP, and includes all SNPs within 0.5*windowsize to the left and right.
    """
    import pandas

    # set up args
    args = kwargs['args']
    mincov = args.mincov
    maxcov = args.maxcov
    windowsize = args.maxcov
    mincount = args.mincount
    pop = kwargs['pop']

    # set up df
    df = snps.loc[chrom,:].copy()
    if isinstance(df, pandas.core.series.Series):
        # if chrom has only one row
        df = pandas.DataFrame(df).T
    locuscol = 'unstitched_locus' if 'unstitched' in kwargs['chromcol'] else 'locus'
    poscol = 'unstitched_pos' if 'unstitched' in kwargs['chromcol'] else 'POS'
    df.index = df[locuscol].tolist()

    windows = {}
    for locus in df.index:
        # locus is the window's ID (the window is centered on the locus)
        windows[locus] = {'measure': args.measure,
                          'snps': {}
                         }

        # find all SNPs within window
        chrom,pos = locus.split('-')
        rows = abs(df[poscol]-int(pos)) <= 0.5*windowsize

        # get info for all SNPs within window
        for snp in df.loc[rows, locuscol]:
            # if snp passes mincount (min count of the minor allele), min/maxcov flags
            # (first part checks AD and RD so I don't have to determine minor allele)
            if all([df[f'{pop}.RD'][snp] >= mincount, df[f'{pop}.AD'][snp] >= mincount,
                    df[f'{pop}.DP'][snp] >= mincov, df[f'{pop}.DP'][snp] <= maxcov]):
                # eucov is the count of all alleles passing quality and mincount filters
                # .DP is the depth coverage for bases >= quality threshold given to varscan (pipeline sets 20)
                # assert mincov <= eucov <= maxcov
                # iscov means eucov passed assertion
                # count_covered = num(iscov) = covercount = num snps in window passing qual AND cov thresholds
                windows[locus]['snps'][snp] = {'eucov': int(df[f"{pop}.DP"][snp]),
                                               'alleles': {'ref': df[f"{pop}.RD"][snp],
                                                           'alt': df[f"{pop}.AD"][snp]
                                                          }}
        windows[locus]['covercount'] = len(windows[locus]['snps'])

    return windows


def write_tmp_file(measures, chrom, pop, statistic, outdir, inputfile, windowsize) -> str:
    """Write window file to /tmp, combine later."""
    import os
    
    bname = os.path.basename(inputfile).replace(".txt", "")

    tmpdir = os.path.join(outdir, 'tmp')
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    file = os.path.join(tmpdir, f"{bname}_{pop}_{statistic}_{windowsize}bp-windows_{chrom}.txt")
    with open(file, 'w') as o:
        o.write(f"chrom\tlocus\tsnpcount\tcoveredFraction\tstatistic ({statistic})\n")
        lines = []
        for locus,dic in measures.items():
            lines.append('\t'.join([chrom, locus, str(dic['covercount']),
                                    str(dic['covfraction']), str(dic[statistic])] ))
        o.write("\n".join(lines))

    return file


def send_windows(*args, **kwargs) -> list:
    """Send each window to varmath.VarianceExactCorrection."""
    import numpy, varmath as vm

    chroms = args[0]

    # gather args, kwargs
    args = kwargs['args']
    vm_kwargs = {'mincov':args.mincov,
                 'maxcov':args.maxcov,
                 'poolsize':kwargs['ploidy'][kwargs['pop']],
                 'mincount':args.mincount,
                 'pidiv_buffer':kwargs['pidiv_buffer'],  # TODO: not needed after deprecating send_get_pidiv_buffer()
                }

    # do the popgen
    files = []
    for chrom in chroms:
        measures = {}
        for locus,window in get_windows(chrom, **kwargs).items():
            # filter by covercount
            measures[locus] = {}
            covfraction = (len(window['snps']) / args.windowsize)
            if covfraction < args.mincovfraction or len(window['snps'])==0:
                # if the fraction of the window occupied by SNPs passing filtered is not high enough
                measures[locus][args.measure] = numpy.nan
            else:
                measures[locus][args.measure] = vm.VarianceExactCorrection(**vm_kwargs).calculate_measure(**window)
            measures[locus]['covercount'] = window['covercount']  # number of snps passing filters in window
            measures[locus]['covfraction'] = window['covercount'] / args.windowsize

        # write the file
        file = write_tmp_file(measures, chrom, kwargs['pop'], args.measure, args.outdir, args.input, args.windowsize)
        files.append(file)

    return files


def send_chrom_to_calculator(lview, **kwargs) -> str:
    """
    Parallelize varmath.VarianceExactCorrection by individual chroms.

    # TODO: make option to parallelize across windows instead of chroms (assuming chroms >> windows)
    """
    import tqdm, math

    # determine how to divy up jobs
    chroms = uni(snps[kwargs['chromcol']])
    jobsize = math.ceil(len(chroms)/1000)

    # send jobs to engines
    numjobs = (len(chroms)/jobsize)+1
    print(ColorText("\nSending %d jobs to engines ..." % numjobs ).bold())
    jobs = []
    chroms_to_send = []
    count = 0
    for chrom in tqdm.tqdm(chroms):
        count += 1
        chroms_to_send.append(chrom)
        if len(chroms_to_send) == jobsize or count == len(chroms):
            jobs.append(lview.apply_async(send_windows, *[chroms_to_send], **kwargs))
#             jobs.append(send_windows(*[chroms_to_send], **kwargs))  # for debugging
            chroms_to_send = []
#     measure_df = pd.concat([pd.read_table(f) for j in jobs for f in j])  # for debugging

    # wait until jobs finish
    watch_async(jobs, phase='send_windows')

    # get tmp file names, read in, concatenate into one df
    print(ColorText("\nGathering results ...").bold())
    measure_df = pd.concat([pd.read_table(f) for j in jobs for f in j.r])

    # write real file
    args = kwargs['args']
    pop = kwargs['pop']
    bname = os.path.basename(args.input).replace(".txt", "")
    # save statistics file
    print(ColorText("\nWriting results to file ...").bold())
    file = op.join(args.outdir, f"{pop}_{args.measure}_{args.windowsize}bp-windows_{bname}.txt")
    measure_df.to_csv(file, sep='\t', index=False)
    # save input arguments
    pkldump(args, os.path.join(args.outdir,
                               f"{pop}_{args.measure}_{args.windowsize}bp-windows_{bname}_ARGS.pkl"))

    return file


def check_pyversion() -> None:
    """Make sure python is 3.6 <= version < 3.8."""
    pyversion = float(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if not pyversion >= 3.6:
        text = f'''FAIL: You are using python {pyversion}. This pipeline was built with python 3.7.
FAIL: use  3.6 <= python version < 3.8
FAIL: exiting pypoolation.py'''
        print(ColorText(text).fail())
        exit()
    if not pyversion < 3.8:
        print(ColorText("FAIL: python 3.8 has issues with the ipyparallel engine returns.").fail())
        print(ColorText("FAIL: use  3.6 <= python version < 3.8").fail())
        print(ColorText("FAIL: exiting pypoolation.py").fail())
        exit()


def main(pidiv_buffer:dict={}):
    from tqdm import tqdm as pbar
    # make sure it's not python3.8
    check_pyversion()

    # get args
    args = get_pars()

    # get ploidy dict
#     global ploidy  # for debugging send_windows()
    ploidy = get_ploidy(args)


    # iterate through pops
    for pop in ploidy:
        # read in VariantsToTable.txt file, filter chroms based on args
        global snps
        snps, chromcol = get_datatable(args, pop)

        # get ipcluster engines
        lview,dview = launch_engines(args.engines, args.profile)

        # attach data, functions, and dict to engines (used in and downstream of send_to_calculate())
        attach_data(snps=snps,
                    ploidy=ploidy,
                    write_tmp_file=write_tmp_file,
                    send_windows=send_windows,
                    send_chrom_to_calculator=send_chrom_to_calculator,
                    get_windows=get_windows,
#                     pidiv_buffer=pidiv_buffer,
                    dview=dview)

        # calculate measure
        file = send_chrom_to_calculator(lview, chromcol=chromcol, args=args, pop=pop,
                                        pidiv_buffer=pidiv_buffer, ploidy=ploidy)
        print(ColorText("\nWrote stats to ").green().__str__() + ColorText(file).bold().green().__str__())
        print(ColorText("\nWrote pypoolation arguments used to ").green().__str__() +
              ColorText(file.replace(".txt", "_ARGS.pkl")).green().bold().__str__())

        # kill ipcluster to avoid mem problems (restart next loop)
        print(ColorText("\nStopping ipcluster ...").bold())
        subprocess.call([shutil.which('ipcluster'), 'stop'])
        
        # remove temporary files
        print(ColorText("\nRemoving temporary files ...").bold())
        tmpdir = os.path.join(args.outdir, 'tmp')
        for f in pbar(os.listdir(tmpdir)):
            os.remove(os.path.join(tmpdir, f))

        print(ColorText("\nDONE!!!").bold().green())

if __name__ == '__main__':
    mytext = ColorText('''
*****************************************************************************
                                CoAdapTree's
        _____         _____             _         _
       |     \\       |     \\           | |     __| |_ O
       |  D  |       |  D  | __    __  | |  __ |__ __|_   __   _ __
       |  __/  _   _ |  __/ /  \\  /  \\ | | /  V| | | | | /  \\ | V  |
       | |    \\ \\/ / | |   | O  || O  || |  C  | | | | || O  || || |
       |_|     \\  /  |_|    \\__/  \\__/ |_| \\_/_| |_| |_| \\__/ |_||_|
               / /             . . 
              /_/           \\        /
                             \\______/
                                 U
                
                Tajima's Pi, Watterson's Theta, Tajima's D
*****************************************************************************''').green().bold().__str__()
    main()

