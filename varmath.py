"""Functions to calculate pi, theata, D."""


import numpy


def get_nbase_buffer(poolsize:int) -> float:
    """
    poolsize - The haploid number of the pool which has been sequenced.

    Returns:
    nbase - the expected number of distinct individuals sequenced
    """

    def _get_pij_matrix(maxcoverage:int, poolsize:int) -> numpy.matrix:
        import numpy, pandas
        jboundary = maxcoverage if maxcoverage < poolsize else poolsize
        matrix = numpy.matrix(pandas.DataFrame(index=range(0, maxcoverage+1, 1),
                                               columns=range(0, jboundary+1, 1)))
        matrix[0,0] = 1
        for i in range(1, maxcoverage+1, 1):
            matrix[i, 0] = 0
        for j in range(1, jboundary+1, 1):
            matrix[0, j] = 0

        for i in range(1, maxcoverage+1, 1):
            for j in range(1, jboundary+1, 1):
                t1 = ((1+poolsize-j)/poolsize)*(matrix[i-1, j-1])
                t2 = (j/poolsize)*(matrix[i-1, j])
                matrix[i,j] = t1 + t2

        return matrix

    mat = _get_pij_matrix(3*poolsize, poolsize)

    # in the perl code, an inequality is used to define minj
    # I can reproduce all numbers when setting minj to poolsize,
        # and using poolsize in the mat look-up
    # minj = poolsize

    nbase = 0
#     for k in range(1, minj+1, 1):v  # <- had I used an inequality to determine minj
    for k in range(1, poolsize+1, 1):
        x = mat[poolsize, k]
        nbase += k * x

    return nbase


def get_an_buffer(n:float, an_buffer={}) -> (float, dict):
    """
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    an_buffer - dict for easy lookup of an(poolsize (->n))
    
    Returns:
    an - value needed for fstar, alphastar, and betastar
    an_buffer - updated dict with key = n, for fast and easy lookup of an
    """

    # TODO: Use a decorator to wrap get_an_buffer() and get_bn_buffer()

    import math

    # in perl, floats were passed to iterators. That outcome is accomplished by .floor
    n = math.floor(n)

    if n in an_buffer.keys():
        return an_buffer[n], an_buffer

    an = 0
    for i in range(1, n, 1):
        an += 1/i

    an_buffer[n] = an

    return an, an_buffer


def get_bn_buffer(n:float, bn_buffer:dict={}) -> float:
    """
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    bn_buffer -  dict for easy lookup of value from get_bn_buffer(poolsize (->n)

    Returns:
    bn -  value needed for betastar
    """

    import math

    # in perl, floats were passed to iterators. That outcome is accomplished by .floor
    n = math.floor(n)

    if n in bn_buffer.keys():  # if I weren't parallelizing, this would help
        return bn_buffer[n]

    bn = 0
    for i in range(1, n, 1):
        bn += 1/(i**2)

    bn_buffer[n] = bn

    return bn


def calculate_fstar(an:float, n:float) -> float:
    """
    an - value needed for fstar, alphastar, and betastar
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    
    Returns:
    fstar - used for alphastar and betastar calculations
    """
    return ((n - 3)/(an*(n - 1) - n))


def get_betastar_calculator(n:float, an_buffer:dict) -> float:
    """
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    an_buffer - dict for easy lookup of an

    Returns:
    bstar -  used in the estimation for the variance of (pi - theta)
    """

    from pypoolation import ColorText

    if not n > 1:
        print(ColorText("\nFAIL: invalid effective coverage; has to be larger than 1").fail())
        print(ColorText("FAIL: exiting pypoolation.py").fail())
        exit()  # TODO: do I care that this should be sys.exit()? If so, should I change everything?

    an, an_buffer = get_an_buffer(n, an_buffer)
    bn = get_bn_buffer(n)
    fs = calculate_fstar(an, n)

    t1 = (fs**2) * (bn - ((2*(n-1)) / ((n-1)**2)))
    st1 = bn * (8/(n-1))
    st2 = an * (4/(n*(n-1)))
    st3 = ((n**3)+12*(n**2)-35*n+18)/(n*((n-1)**2))
    t2 = fs*(st1-st2-st3)
    t3 = bn * (16/(n*(n-1)))
    t4 = an * (8/((n**2)*(n-1)))
    st4 = 2*(n**4+ 110*(n**2)-255*n+126)
    st5= 9*(n**2)*((n-1)**2)
    t5 = st4/st5

    bstar = (t1 + t2 - t3 + t4 + t5)

    return bstar


def get_alphastar_calculator(n:float) -> (float, dict):
    """
    poolsize (->n) - The haploid number of the pool which has been sequenced.

    Returns:
    astar -  used in the estimation for the variance of (pi - theta)
    an_buffer - updated dict with key = n, for fast and easy lookup of an
    """
    
    an, an_buffer = get_an_buffer(n)
    fs = calculate_fstar(an, n)
    
    t1 = (fs**2)*(an-(n/(n-1)))
    st1 = an * ( (4*(n+1)) / ((n-1)**2) )
    st2 = 2 * ((n+3)/(n-1))
    t2 = fs * (st1-st2)
    t3 = an * ( (8*(n+1))/(n*((n-1)**2)) )
    t4 = ((n**2)+n+60)/(3*n*(n-1))

    astar = (t1 + t2 - t3 + t4)

    return astar, an_buffer


def get_ddivisor(n:int, mincoverage:int, snps:dict, theta:int, nbase_buffer:dict={}) -> float:
    """
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    mincoverage - The minimum coverage of a site, input by user. Sites with a lower coverage will not be considered.
    snps - dict from get_windows(). keys = locus names within a window. value = dict, with keys:
        eucov (->M) - value = total depth of coverage
        alleles - dict, with keys = ['ref', 'alt'], values = {'REF_depth', 'ALT_depth'}
    theta - result from get_theta_calculator() - Watterson's theta for the window
    nbase_buffer - dict with key=n for fast and easy lookup of nbase - the expected number of distince individuals sequenced

    Returns:
    div - denominator of Tajima's D - estimator of standard deviation of (Tajima's pi - Tajima's theta)
    """

    import numpy

    snpcount = len(snps.keys())
    averagen = get_nbase_buffer(n)
#     # This's how it appears in popoolation code, is betastar_calc right for returning alphastar??
#     alphastar = get_betastar_calculator(averagen)
    # since it ^ doesn't make sense I'll use get_alphastar_calculator():
    alphastar, an_buffer = get_alphastar_calculator(averagen)
    betastar = get_betastar_calculator(averagen, an_buffer)
    div = numpy.sqrt(alphastar/snpcount)*theta + (betastar*(theta**2))

    return div


def get_thetadiv_buffer(b:int, n:int, M:int, thetadiv_buffer:dict, amnm_buffer:dict) -> (float, dict, dict):
    """
    b - The minimum count of the minor allele, input by user.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    eucov (->M) - value = total depth of coverage
    thetadiv_buffer - dict with key b:n:M permutation, for fast and easy lookup ...
        of value from get_thetadiv_buffer(b, poolsize (->n), eucov (->M), thetadiv_buffer, amnm_buffer)
    amnm_buffer - dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()

    Returns:
    div -  value used to calculate it's reciprical in the summation of thetasum
    thetadiv_buffer - updated dict with key b:n:M permutation, for fast and easy lookup ...
        of value from get_thetadiv_buffer(b, poolsize (->n), eucov (->M), thetadiv_buffer, amnm_buffer)
    amnm_buffer - updated dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    """

    key = f"{b}-{n}-{M}"

    if key in thetadiv_buffer:
        return thetadiv_buffer[key], thetadiv_buffer, amnm_buffer

    div = 0
    for m in range(b, M+1-b, 1):
        buff, amnm_buffer = get_aMnm_buffer(M, n, m, amnm_buffer)
        div += buff
    thetadiv_buffer[key] = div

    return div, thetadiv_buffer, amnm_buffer


def get_theta_calculator(b:int, n:int, snps:dict, amnm_buffer:dict, thetadiv_buffer:dict={}) -> float:
    """Calculate Watterson's theta statistic.
    
    b - The minimum count of the minor allele, input by user.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    snps - dict from get_windows(). keys = locus names within a window. value = dict, with keys:
        eucov (->M) - value = total depth of coverage
        alleles - dict, with keys = ['ref', 'alt'], values = {'REF_depth', 'ALT_depth'}
    amnm_buffer - dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    thetadiv_buffer - dict with key b:n:M permutation, for fast and easy lookup ...
        of value from get_thetadiv_buffer

    Returns:
    thetasum - estimator of Watterson's Theta
    """

    thetasum = 0
    for snp,snpdict in snps.items():
        div, thetadiv_buffer, amnm_buffer = get_thetadiv_buffer(b, n, snpdict['eucov'], thetadiv_buffer, amnm_buffer)
        thetasum += (1 / div)

    return thetasum


def binomial_term(M:int, n:int, m:int, k:int) -> float:
    """
    eucov (->M) - value = total depth of coverage
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    m - iterator from get_pidiv_buffer(), for m in range(b, M+1-b, 1)
    k - iterator from get_aMnm_buffer(), for k in range(1, n, 1)

    Returns:
    # The probability of having a first allele count of m in M reads from ...
        a pool of n with first allele count of k (m is the allele count ...
        in the reads, k is the allele count in the pool)
    """

    from scipy.special import comb

    try:
        assert M > m
    except AssertionError as e:
        print ('M = %s, k = %s' % (M,k))
        raise AssertionError

    val = comb(M, m)
    t1 = (k/n)**m
    t2 = ((n-k)/n)**(M-m)

    return val*t1*t2


def get_aMnm_buffer(M:int, n:int, m:int, amnm_buffer:dict) -> (float, dict):
    """
    eucov (->M) - value = total depth of coverage
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    m - iterator from get_pidiv_buffer(), for m in range(b, M+1-b, 1)

    Returns:
    toret - value from amnm_buffer[key]
    amnm_buffer - updated dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    """

    key = f"{M}-{n}-{m}"

    if key in amnm_buffer:
        return  amnm_buffer[key], amnm_buffer

    toret = 0
    for k in range(1, n, 1):
        t1 = binomial_term(M, n, m, k)
        t1 *= 1/k
        toret += t1

    amnm_buffer[key] = toret
    
    return toret, amnm_buffer


def get_pidiv_buffer(b:int, n:int, M:int, amnm_buffer:dict, pidiv_buffer:dict={}) -> (float, dict, dict):
    """
    b - The minimum count of the minor allele, input by user.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    eucov (->M) - value = total depth of coverage
    amnm_buffer - dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    pidiv_buffer - dict with key b:n:M permutations, for fast and easy lookup ...
        of value from get_pidiv_buffer()

    Returns:
    div - value from pidiv_buffer[M:n:m permutation]
    amnm_buffer - updated dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    pidiv_buffer - updated dict with key b:n:M permutations, for fast and easy lookup ...
        of value from get_pidiv_buffer()
    """

    key = f"{b}-{n}-{M}"

    if key in pidiv_buffer:
        return pidiv_buffer[key], pidiv_buffer, amnm_buffer

    div = 0
    amnm_buffer = {}
    for m in range(b, M+1-b, 1):
        term1 = (2*m*(M-m))/(M*(M-1))
        buff, amnm_buffer = get_aMnm_buffer(M, n, m, amnm_buffer)
        term1 *= buff
        div += term1

    pidiv_buffer[key] = div

    return div, pidiv_buffer, amnm_buffer


def get_pi_calculator(b:int, n:int, snps:dict, pidiv_buffer:dict={}, amnm_buffer:dict={}) -> (float, dict):
    """Calculate Tajima's pi.
    
    b - The minimum count of the minor allele, input by user.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    snps - dict from get_windows(). keys = locus names within a window. value = dict, with keys:
        eucov (->M) - value = total depth of coverage
        alleles - dict, with keys = ['ref', 'alt'], values = {'REF_depth', 'ALT_depth'}
    pidiv_buffer - dict with key b:n:M permutations, for fast and easy lookup ...
        of value from get_pidiv_buffer()
    amnm_buffer - dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()

    Returns:
    pisum - Tajima's pi for SNPs in window
    amnm_buffer - updated dict with value of M:n:m permutation, for fast and easy lookup ...
        of value from get_aMnm_buffer()
    """

    pisum = 0
    for snp,snpdict in snps.items():
        M = snpdict['eucov']

        pi_snp = 1
        for allele,count in snpdict['alleles'].items():
            pi_snp -= (count/M)**2

        pi_snp *= M/(M-1)

        div, pidiv_buffer, amnm_buffer = get_pidiv_buffer(b, n, M, amnm_buffer, pidiv_buffer)

        if div == 0:  # can't div by zero
            print('bad 0 ', b, n, M)

        pi_snp /= div

        pisum += pi_snp

    return pisum, amnm_buffer


def get_D_calculator(b:int, n:int, mincoverage:int, snps:dict, pidiv_buffer:dict={}) -> (float, dict):
    """Calculate Tajima's D (corrected).
    
    b - The minimum count of the minor allele.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    mincoverage - The minimum coverage of a site, input by user. Sites with a lower coverage will not be considered.
    snps - dict from get_windows(). keys = locus names within a window. value = dict, with keys:
        eucov (->M) - value = total depth of coverage
        alleles - dict, with keys = ['ref', 'alt'], values = {'REF_depth', 'ALT_depth'}
    pidiv_buffer - dict with key b:n:M permutations, for fast and easy lookup ...
        of value from get_pidiv_buffer()
    
    # notes:
    - passing pidiv_buffer from get_D_calculator was from when I
        parallelized pidiv_buffer before measure calc. Without this,
        there is no need to have it passed into get_D_calculator

    Returns:
    d - Tajima's D estimate for window
    """

    pi, amnm_buffer = get_pi_calculator(b, n, snps, pidiv_buffer)  # pidiv_buffer, amnm_buffer initiated
    theta = get_theta_calculator(b, n, snps, amnm_buffer)  # thetadiv_buffer uses amnm_buffer
    ddivisor = get_ddivisor(n, mincoverage, snps, theta)

    if pi - theta == 0:
        return 0

    d = (pi - theta) / ddivisor

    return d


class VarianceExactCorrection:
    """Calculate popgen statistic from window.

    mincoverage - The minimum coverage of a site. Sites with a lower coverage will not be considered.
    maxcoverage - The maximum coverage of a site. Sites with greater coverage will not be considered.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
    mincount (->b) - The minimum count of the minor allele for a SNP to be considered.
    pidiv_buffer - dict with key b:n:M permutations, for fast and easy lookup ...
        of value from get_pidiv_buffer()
    
    Returns:
    VarianceExactCorrection class object
    """

    from pypoolation import ColorText
    
    def __init__(self, **kwargs):
        self.mincoverage = kwargs['mincov']
        self.maxcoverage = kwargs['maxcov']
        self.n = kwargs['poolsize']
        self.b = kwargs['mincount']
        self.pi = get_pi_calculator
        self.theta = get_theta_calculator
        self.d = get_D_calculator
        self.pidiv_buffer = kwargs['pidiv_buffer']
#         global pidiv_buffer
#         pidiv_buffer = kwargs['pidiv_buffer']


    def calculate_measure(self, **kwargs):
        """Calculate and return specific measure (pi, theta, D) for a window."""
        measure = kwargs['measure'].lower()
        if measure == 'pi':
            return self._calculate_pi(snps=kwargs['snps'], covercount=kwargs['covercount'])
        elif measure == 'theta':
            return self._calculate_theta(snps=kwargs['snps'], covercount=kwargs['covercount'])
        elif measure == 'd':
            return self._calculate_d(snps=kwargs['snps'], covercount=kwargs['covercount'],
                                     mincount=self.b, pdiv_buffer=self.pidiv_buffer)
        else:
            print(ColorText("FAIL: You have not selected a valid option for popgen statistic.").fail())
            print(ColorText("FAIL. The options are: pi, theta, D (case insensitive).").fail())
            print(ColorText("exiting pypoolation.py").fail())
            exit()

    def _calculate_pi(self, **kwargs):
        """Calculate and return the average pi over the window.
        Average is only used when calculating pi in isolation.
        """
        if kwargs['covercount'] == 0:
            return 0
        pi_sum, amnm_buffer = self.pi(b=self.b, n=self.n, snps=kwargs['snps'], pidiv_buffer=self.pidiv_buffer)

        return pi_sum / kwargs['covercount']

    def _calculate_theta(self, **kwargs):
        """Calculate and return the average theta over the window.
        Average is only used when calculating pi in isolation.
        """
        if kwargs['covercount'] ==0:
            return 0
        theta = self.theta(b=self.b, n=self.n, snps=kwargs['snps'], amnm_buffer={})

        return theta/kwargs['covercount']

    def _calculate_d(self, **kwargs):
        """Calculate and return Tajima's D for a window."""
        # TODO: from popoolation: _calculate_d : die if 3*$mincoverage < $poolsize
            # die "Corrected Tajima's D error\n
            # Poolsize >> mincoverage (as internal aproximation: 3 * minimumcoverage < poolsize)" unless 3*$mincoverage < $poolsize;
        if kwargs['covercount'] == 0:
            return 0
        d = self.d(b=self.b, n=self.n, mincoverage=self.mincoverage,
                   snps=kwargs['snps'], pidiv_buffer=self.pidiv_buffer)

        return d
