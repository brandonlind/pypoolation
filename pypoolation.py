"""A python version of some scripts from popoolation.

TODO: whichpops default
TODO: auto stop ipcluster if __main__ dies
TODO: split up main caller (pypoolation) with misc function
TODO: ploidy pkl has first key of pool_name - how do I figure that out computationally? - ask the user for input
TODO: remove ploidy key 'JP_pooled'
TODO: do I assert pop arg is in ploidy[pool_chosen_by_user]?
TODO: remove all "global [object]"

"""

import os, sys, argparse, shutil, subprocess, pandas as pd, threading, ipyparallel, time, signal
import itertools, pickle, numpy as np, math
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
    """Load object from a .pkl file"""
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

def make_jobs(inputs:list, cmd, lview) -> list:
    """Send each arg from inputs to a function command; async."""
    print(f"making jobs for {cmd.__name__}")
    jobs = []
    for arg in tnb(inputs):
        jobs.append(lview.apply_async(cmd, arg))
    return jobs


def watch_async(jobs:list, phase=None) -> None:
    """Wait until jobs are done executing, show progress bar."""
    from tqdm import trange
    print(ColorText(f"\nWatching {len(jobs)} {phase} jobs ...").bold())
    for i in trange(len(jobs)):
        count = 0
        while count < (i+1):
            count = 0
            for j in jobs:
                if j.ready():
                    count += 1


class ColorText():
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
    """Parse input flags."""
    parser = argparse.ArgumentParser(description=mytext,
                                     add_help=True,
                                     formatter_class=argparse.RawTextHelpFormatter)
    requiredNAMED = parser.add_argument_group('required arguments')
    requiredNAMED.add_argument("-i", "--input",
                               required=True,
                               default=None,
                               dest="input",
                               type=str,
                               help="/path/to/VariantsToTable_output.txt")
    requiredNAMED.add_argument("-o","--outdir",
                               required=True,
                               default=None,
                               dest="outdir",
                               type=str,
                               help="/path/to/pypoolation_output.txt")
    requiredNAMED.add_argument("-m","--measure",
                               required=True,
                               default=None,
                               dest="measure",
                               type=str,
                               help='''pi OR theta OR D. If D is selected, both pi and theta will
also be output.''')
    requiredNAMED.add_argument("-p","--ploidy",
                               required=True,
                               default=None,
                               dest="ploidyfile",
                               type=str,
                               help="/path/to/the/ploidy.pkl file output by the VarScan pipeline.")
    requiredNAMED.add_argument("-e","--engines",
                               required=True,
                               default=None,
                               dest="engines",
                               type=int,
                               help="The number of multiprocessing engines that will be launched.")
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
and max-coverage (fraction = num snps / windowsize). Default=0.0''')
    parser.add_argument("--window-size",
                        required=False,
                        default=1000,
                        dest="windowsize",
                        type=str,
                        help='''The size of the sliding window. Windows are centered on
individual SNPs, so 0.5*windowsize bps will be on either
side of a SNP. Default=1000bp''')
    parser.add_argument("--step-size",
                        required=False,
                        default=1000,
                        dest="stepsize",
                        type=int,
                        help='''size of one sliding window step. If this number is equal
to the --window-size the sliding window will be non 
overlapping (jumping window). Default 1000bp''')
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
        print(ColorText(f"FAIL: the directory for the output file(s) does not exist.").bold().fail())
        print(ColorText(f"FAIL: please create this directory: %s" % args.outdir).bold().fail())
        print(ColorText("exiting pypoolation.py").bold().fail())  # TODO: check to make sure I'm consistent with .bold().fail()
        exit()
    nopath = []
    for x in [args.input, args.ploidyfile]:
        if not op.exists(x):
            nopath.append(x)
    if len(nopath) > 0:
        print(ColorText("FAIL: The following path%s do not exist:" % "s" if len(nopath) > 1 else "").fail())
        [print(ColorText("FAIL: %s" % f).fail()) for f in nopath]
    
    print('args = ', args)
    return args


def askforinput(msg='Do you want to proceed?', tab='', newline='\n'):
    "Ask for input; if msg is default and input is no, exit."
    while True:
        inp = input(ColorText(f"{newline}{tab}INPUT NEEDED: {msg} \n{tab}(yes | no): ").warn().__str__()).lower()
        if inp in ['yes', 'no']:
            if inp == 'no' and msg=='Do you want to proceed?':
                print(ColorText('exiting %s' % sys.argv[0]).bold().fail())
                exit()
            break
        else:
            print(ColorText("Please respond with 'yes' or 'no'").fail().bold())
    return inp


def get_datatable(args, pop, chroms=None):
    """
    Load --input datatable.txt.
    """
    print(ColorText(f"\nReading in SNP datatable for {pop}...").bold())
    # read in path
    df = pd.read_table(args.input)

    # see if user supplied specific chromosomes/contigs to analyze
    if args.chromfile is not None:
        if op.exists(args.chromfile):
            chroms = read(args.chromfile)
        else:
            print(ColorText(f"\nFAIL: chromsome file does not exist: {args.chromfile}").bold().fail())
            print(ColorText("exiting pypoolation.py").bold().fail())
            exit()

    # determine which column to use for parallelization
    chromcol = 'unstitched_chrom' if 'unstitched_chrom' in df.columns else 'CHROM'
    if chromcol not in df.columns:
        print(ColorText(f"FAIL: {chromcol} is not in the input columns").bold().fail())
        print(ColorText("exiting pypoolation.py").fail.bold())
        exit()
    print(f'\tUsing the {chromcol} column to parallelize data')
    
    # reduce df to pop that matters
    keepcols = [col for col in df.columns if '.' not in col or pop in col]
    df = df[keepcols].copy()
    
    # reduce rows
    df = df[df[f'{pop}.FREQ']==df[f'{pop}.FREQ']].copy()  # only keep filtered SNPs (filtered across all pops)
    df = df[df[f'{pop}.FREQ']!= '0%'].copy()  # remove SNPs that are invariable in the current pop
    
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
    count = 0
    while len(lview) != engines:
        if count % 20 == 0:
            # if waiting too long..
            print(ColorText("\tFAIL: waited too long for engines. Exiting pypoolation.").fail())
            print(ColorText("\tFAIL: Make sure that if any cluster is running, the -e arg matches the number of engines.").fail())
            print(ColorText("\tFAIL: In some cases, not all expected engines can start on a busy server.").fail())
            print(ColorText("\tFAIL: pypoolation found %s engines, with -e set to %s" % (len(lview), engines)).fail())
            answer = askforinput(msg='Would you like to continue with %s engines? (choosing no will wait another 30 seconds)' % len(lview), tab='\t', newline='')
            if answer == 'yes':
                break
        try:
            lview,dview = get_client(profile=profile)
        except (OSError, ipyparallel.error.NoEnginesRegistered, ipyparallel.error.TimeoutError) as e:
            lview = []
            pass
        time.sleep(2)
        count += 1
    return lview,dview


def launch_engines(engines, profile):
    """Launch ipcluster with engines under profile."""
    print(ColorText(f"\nLaunching ipcluster with {engines} engines...").bold())
    
    def _launch(engines, profile):
        subprocess.call([shutil.which('ipcluster'), 'start', '-n', str(engines), '--daemonize'])
    
    # first see if a cluster has already been started
    started = False
    try:
        print("\tLooking for exising engines ...")
        lview,dview = get_client(profile=profile)
        if len(lview) != engines:
            lview,dview = wait_for_engines(engines, profile)
        started = True
    except (OSError, ipyparallel.error.NoEnginesRegistered, ipyparallel.error.TimeoutError) as e:
        print("\tNo engines found ...")
        pass
    
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
    print(ColorText("\nAdding data to engines ...").bold())
    print(ColorText("\tWARN: Watch available mem in another terminal window: 'watch free -h'").warn())
    print(ColorText("\tWARN: If available mem gets too low, kill engines and restart pypoolation with fewer engines: 'ipcluster stop'").warn())
    for key,value in kwargs.items():
        if key != 'dview':
            print(f'\tLoading {key} ({value.__class__.__name__}) to engines')
            kwargs['dview'][key] = value
    return None
        

def uni(lst):
    """Return unqiue items from lst."""
    return list(set(lst))

def get_combos(pop, df, ploidy, dview):
    """Get unique combinations of minor_count (b), pop ploidy (n), and coverage (M)."""
    # get unique mincounts (b)  # TODO: don't worry about anthing less than args.mincount
    df[f'{pop}.minor_count'] = np.nan
    minors = df[f'{pop}.AD'] < df[f'{pop}.RD']
    df.loc[minors, f'{pop}.minor_count'] = df.loc[minors, f'{pop}.AD']
    df.loc[~minors, f'{pop}.minor_count'] = df.loc[~minors, f'{pop}.RD']
    uni_b = uni([x for x in df[f'{pop}.minor_count'].astype(int) if x==x and x!=0])
    print('uni_b = ', uni_b)

    # get unique poolsizes (from ploidy)
    uni_n = uni(ploidy.values())

    # get unique coverages
    uni_M = uni([x for x in df[f'{pop}.DP'].astype(int) if x==x])

    # get unique combinations
    uni_combos = list(itertools.product(range(len(uni_b)),
                                        range(len(uni_M)),
                                        range(len(uni_n))))
    
    # attache data to engines
    attach_data(uni_b=uni_b, uni_M=uni_M, uni_n=uni_n, dview=dview)
    
    return uni_combos, uni_b, uni_M, uni_n


def pidiv_iterator(chunk, uni_b, uni_M, uni_n):
    """For each chunk (a unique combo of b, M, and n), calculate pidiv_buffer."""
    import varmath
    
    pi_buffer = {}
    for b_idx, M_idx, n_idx in chunk:
        b = int(uni_b[int(b_idx)])
        M = int(uni_M[int(M_idx)])
        n = int(uni_n[int(n_idx)])
        key = f"{b}-{n}-{M}"
        div, buffdict = varmath.get_pidiv_buffer(b, n, M)
        pi_buffer[key] = div
    return pi_buffer


def send_get_pidiv_buffer(**kwargs):
    """Parallelize varmath.get_pidiv_buffer()."""
    
    # get unique combos of b, n, and M
    uni_combos, uni_b, uni_M, uni_n = get_combos(kwargs['pop'], kwargs['snps'], kwargs['ploidy'], kwargs['dview'])
    for name,lst in zip(['uni_combos', 'uni_b', 'uni_M', 'uni_n'],[uni_combos, uni_b, uni_M, uni_n]):
        print(name, len(lst))
    print(uni_combos[:5])
    
    # send out jobs
    maxjobs = math.ceil(len(uni_combos)/kwargs['engines'])
    print('maxjobs = ', maxjobs)
    jobs = []
    chunk = []
    for i,combo in enumerate(uni_combos):
        chunk.append(combo)
        if len(chunk) == maxjobs or (i+1) == len(uni_combos):
            jobs.append(kwargs['lview'].apply_async(pidiv_iterator, *(chunk, uni_b, uni_M, uni_n)))
#             jobs.append(pidiv_iterator(chunk, uni_b, uni_M, uni_n))
            chunk = []

    # wait for jobs to finish
    watch_async(jobs, 'pidiv_iterator')
    
    # get returns
    buffdict = {}
    for j in jobs:
        print('len(j.r) = ', len(j.r))
        print('type(j.r) = ', type(j.r))
        buffdict.update(j.r)
    print('buffdict = ', list(buffdict.keys())[:5])

    return buffdict
    


def get_ploidy(ploidyfile, pops):
    """Get the ploidy of the populations of interest, reduce ploidy pkl."""
    global ploidy
    ploidy = pklload(ploidyfile)['JP_pooled']  # TODO: remove this
     
    for pop in list(ploidy.keys()):
        if not pop in pops:
            ploidy.pop(pop)
    return ploidy


def get_windows(chrom, **kwargs):
    """Get info for all SNPs meeting criteria for a given window size on a specific chrom.
    
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


def write_tmp_file(measures, chrom, pop, statistic):
    """Write window file to /tmp, combine later.
    
    # TODO: include flag values in outfile name
    """
    import tempfile, os

    file = os.path.join(tempfile.gettempdir(), f"{chrom}_{statistic}_{pop}.txt")
    with open(file, 'w') as o:
        o.write(f"chrom\tlocus\tsnpcount\tcoveredFraction\tstatistitc ({statistic})\n")
        lines = []
        print ("chrom = ", chrom)
        for locus,dic in measures.items():
            lines.append('\t'.join([chrom, locus, str(dic['covercount']),
                                    str(dic['covfraction']), str(dic[statistic])] ))
        o.write("\n".join(lines))
        
    return file


def send_windows(*args, **kwargs):
    """
    Send each window to varmath.VarianceExactCorrection.
    
    # TODO : make sure this will work for other measures other than D
    """
    import numpy, varmath as vm
    chrom = args[0]
    
    # gather kwargs
    args = kwargs['args']
    print('lucky args = ', args)
    vm_kwargs = {'mincov':args.mincov,
                 'maxcov':args.maxcov,
                 'poolsize':ploidy[kwargs['pop']],
                 'mincount':args.mincount,
                 'pidiv_buffer':kwargs['pidiv_buffer'],  # TODO: it's already in globals, can I just use that instead of passing?
                }

    # do the popgen
    measures = {}  # TODO: I don't need to return a dict, I just need to write to a file.
    for locus,window in get_windows(chrom, **kwargs).items():
        # filter by covercount
        measures[locus] = {}
        if (len(window) / args.windowsize) < args.mincovfraction:
            # if the fraction of the window occupied by SNPs passing filtered is not high enough
            measures[locus][args.measure] = numpy.nan
        else:
            print('window keys = ', window.keys())
            print('window["snps"].keys() = ', window["snps"].keys())
            measures[locus][args.measure] = vm.VarianceExactCorrection(**vm_kwargs).calculate_measure(**window)
        measures[locus]['covercount'] = window['covercount']  # number of snps passing filters in window
        measures[locus]['covfraction'] = window['covercount'] / args.windowsize
    
    # write the file
    file = write_tmp_file(measures, chrom, kwargs['pop'], args.measure)
            
    return file


def send_chrom_to_calculator(snps, lview, **kwargs):
    """
    Parallelize varmath.VarianceExactCorrection by individual chroms.
    
    # TODO: make option to parallelize across windows instead of chroms (assuming chroms >> windows)
    """
    jobs = []
    for chrom in uni(snps[kwargs['chromcol']]):
        jobs.append(lview.apply_async(send_windows, *[chrom], **kwargs))
#         jobs.append(send_windows(*[chrom], **kwargs))

    # wait until jobs finish
    watch_async(jobs, phase='send_windows')
    
    # get tmp file names, read in, concatenate into one df
    measure_df = pd.concat([pd.read_table(j.r) for j in jobs])
    
    # write real file
    # TODO: incorporate pool name once user defines for ploidy.pkl
    args = kwargs['args']
    pop = kwargs['pop']
    file = op.join(args.outdir, f"{args.measure}_{pop}.txt")
    measure_df.to_csv(file, sep='\t', index=False)
    # TODO: delete tmp files?
    return file
    

def main():
    # TODO: check that current dir is in the pythonpath so no import errors
#     check_pythonpath()

    # get args
    args = get_pars()
    
    # get ploidy dict
    ploidy = get_ploidy(args.ploidyfile, args.whichpops)
    
    
    # iterate through pops
    for pop in ploidy:
        # read in VariantsToTable.txt file, filter chroms based on args
        global snps
        snps, chromcol = get_datatable(args, pop)

        # get ipcluster engines
        lview,dview = launch_engines(args.engines, args.profile)

        # attach data on all engines
        attach_data(snps=snps, ploidy=ploidy, dview=dview)

        # get pidiv_buffer, attach data to all engines
#         global pidiv_buffer
        pidiv_buffer = send_get_pidiv_buffer(snps=snps, lview=lview, chromcol=chromcol,
                                             ploidy=ploidy, pop=pop, engines=args.engines,
                                             dview=dview)
        print('\tlen(pidiv_buffer) = ', len(pidiv_buffer))
        
        # attach functions and dict to engines (used in and downstream of send_to_calculate())
        attach_data(write_tmp_file=write_tmp_file,
                    send_windows=send_windows,
                    send_chrom_to_calculator=send_chrom_to_calculator,
                    get_windows=get_windows,
                    pidiv_buffer=pidiv_buffer,
                    dview=dview)

        # calculate measure
        file = send_chrom_to_calculator(snps, lview,
                                        chromcol=chromcol, args=args, pop=pop,
                                        pidiv_buffer=pidiv_buffer)
        print(ColorText("\nWrote stats to %s" % file).bold())

        # kill ipcluster to avoid mem problems (restart next loop)
        print(ColorText("\n\tStopping ipcluster ...").bold())
        subprocess.call([shutil.which('ipcluster'), 'stop'])

    

    


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
               / /              . . 
              /_/            \\        /
                              \\______/
                                  U
                
                Tajima's Pi, Watterson's Theta, Tajima's D
*****************************************************************************''').green().bold().__str__()
    main()

