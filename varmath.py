"""Functions to calculate pi, theata, D."""
import numpy, typing


def get_nbase_buffer(poolsize:int, cov:int) -> float:
    def _get_pij_matrix(maxcoverage:int, poolsize:int) -> numpy.matrix:
        import numpy
        jboundary = maxcoverage if maxcoverage < poolsize else poolsize
        matrix = numpy.matrix(pd.DataFrame(index=range(0, maxcoverage+1, 1),
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

    # this is the long version of the one-line if/else statment immediately afterwards
#     if poolsize >= cov:
#         minj = cov if cov > poolsize else poolsize
#     else:
#         minj = cov if cov <= poolsize else poolsize

    # this inequality is different than perl version
        # but with mat[poolsize, k] below (instead of mat[cov][k] in perl) it works
    minj = poolsize if poolsize != cov else cov

    nbase = 0
    for k in range(1, minj+1, 1):
        x = mat[poolsize, k]
        nbase += k * x

    return nbase


def get_an_buffer(n:int):


def get_betastar_calculator(n:int, ):
    from pypoolation import ColorText

    if not n > 1:
        print(ColorText("\nFAIL: invalid effective coverage; has to be larger than 1").fail())
        print(ColorText("FAIL: exiting pypoolation.py").fail())  # TODO: did I check bold convention of ColorText().fail()?
        exit()  # TODO: do I care that this should be sys.exit()? If so, should I change everything?

    if n in bstar_buffer.keys():
        return bstar_buffer[n]

    an = get_an_buffer(n)  # TODO: write fxn()
    bn = get_bn_buffer(n)  # TODO: write fxn()
    myfs = calculate_fstar(an, n)

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
    
    bstar_buffer[n] = bstar
        
    return bstar


def get_ddivisor(n:int, mincoverage:int, snps:dict, theta:int, nbase_buffer:dict={}) -> (float, dict):
    snpcount = len(snps.keys())
    averagen = get_nbase_buffer(n, mincoverage)
    alphastar = get_betastar_calculator(averagen)  # TODO: This's how it appears in popoolation code, is betastar_calc right??
    betastar = get_betastar_calculator(averagen)
    div = (alphastar/snpcount)*theta + (betastar*(theta**2))
    pass


def get_thetadiv_buffer(b:int, n:int, M:int, thetadiv_buffer:dict) -> (float, dict):
    key = f"{b}-{n}-{M}"
    if key in thetadiv_buffer:
        return thetadiv_buffer[key], thetadiv_buffer
    div = 0
    for m in range(b, M+1-b, 1):
        div += get_aMnm_buffer(M, n, m)
    thetadiv_buffer[key] = div
    return div, thetadiv_buffer


def get_theta_calculator(b:int, n:int, snps:dict, thetadiv_buffer:dict={}) -> float:
    thetasum = 0
    for snp,snpdict in snps.items():
        div, thetadiv_buffer = get_thetadiv_buffer(b, n, snpdict['eucov'], thetadiv_buffer)
        thetasum += (1 / div)
    return thetasum


def binomial_term(M:int, n:int, m:int, k:int):
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


def get_aMnm_buffer(M:int, n:int, m:int, amnm_buffer={}) -> float:
    # TODO: make sure this dict is being passed back and forth
    key = f"{M}-{n}-{m}"
    if key in amnm_buffer:
        return  amnm_buffer[key]
    toret = 0
    for k in range(1, n, 1):
        t1 = binomial_term(M, n, m, k)
        t1 *= 1/k
        toret += t1
    amnm_buffer[key] = toret
    return toret


def get_pidiv_buffer(b:int, n:int, M:int, pidiv_buffer:dict={}) -> (float, dict):
    """
    b - The minimum count of the minor allele.
    n - poolsize
    M - coverage
    """
    key = f"{b}-{n}-{M}"
    if key in pidiv_buffer:
        return pidiv_buffer[key], pidiv_buffer
    print('weird key = ', key)
    div = 0
    amnm_buffer = {}  # TODO: is this a good spot for this?
    for m in range(b, M+1-b, 1):
        term1 = (2*m*(M-m))/(M*(M-1))
        buff = get_aMnm_buffer(M, n, m)
        term1 *= buff
        div += term1
    pidiv_buffer[key] = div
    return div, pidiv_buffer


def get_pi_calculator(b:int, n:int, snps:dict, pidiv_buffer:dict={}) -> float:
    """Calculate Tajima's pi.
    
    b - The minimum count of the minor allele.
    n - poolsize
    snps - dict with key=snp, value=dict with keys of coverage and allele counts
    """
    pisum = 0
    for snp,snpdict in snps.items():
        print('keys = ', snps.keys())
        M = snpdict['eucov']
        pi_snp = 1
        for allele,count in snpdict['alleles'].items():
            pi_snp -= (count/M)**2
        pi_snp *= M/(M-1)
        div, pidiv_buffer = get_pidiv_buffer(b, n, M, pidiv_buffer)
        print('type(div) = ', type(div), ' type(pidiv_buffer) = ', type(pidiv_buffer))
        print('div = ', div)
        if div == 0:
            print('bad 0 ', b, n, M)
        pi_snp /= div
        
        pisum += pi_snp

    return pisum


def get_D_calculator(b:int, n:int, mincoverage:int, snps:dict, pidiv_buffer:dict={}) -> (float, dict):
    """Calculate Tajima's D (corrected).
    
    b = The minimum count of the minor allele.
    n = poolsize
    """
    pi = get_pi_calculator(b, n, snps, pidiv_buffer)
    theta = get_theta_calculator(b, n, snps)
    ddivsor = get_ddivisor(n, mincoverage, snps, theta)
    
    if pi - theta == 0:
        return 0
    d = (pi - theta) / ddivisor
    return d, pidiv_buffer

class VarianceExactCorrection:
    """Calculate popgen statistic from window.
    
    mincount (->b) - The minimum count of the minor allele.
    mincoverage - The minimum coverage of a site. Sites with a lower coverage will not be considered.
    maxcoverage - The maximum coverage of a site. Sites with greater coverage will not be considered.
    poolsize (->n) - The haploid number of the pool which has been sequenced.
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
        measure = kwargs['measure'].lower()
        if measure == 'pi':
            return self._calculate_pi(snps=kwargs['snps'], covercount=kwargs['covercount'])
        elif measure == 'theta':
            return self._calculate_theta(snps=kwargs['snps'], covercount=kwargs['covercount'])
        elif measure == 'd':
            return self._calculate_d(snps=kwargs['snps'], covercount=kwargs['covercount'],
                                     mincount=self.b, pdiv_buffer=self.pidiv_buffer)
        else:
            print(ColorText('''FAIL: You have not selected a valid option for popgen statistic.
FAIL. The options are: pi, theta, D (case insensitive).
exiting pypoolation.py''').fail().bold())
            exit()

    def _calculate_pi(self, **kwargs):
        if kwargs['covercount'] == 0:
            return 0
        pi_sum = self.pi(b=self.b, n=self.n, snps=kwargs['snps'], pidiv_buffer=self.pidiv_buffer)
        
        return pi_sum / kwargs['covercount']

    def _calculate_theta(self, **kwargs):
        if kwargs['covercount'] ==0:
            return 0
        theta = self.theta(b=self.b, n=self.n, snps=kwargs['snps'])
        
        return theta/kwargs['covercount']

    def _calculate_d(self, **kwargs):
        # TODO: from popoolation: _calculate_d : die if 3*$mincoverage < $poolsize
        if kwargs['covercount'] == 0:
            return 0
        d = self.d(b=self.b, n=self.n, mincoverage=self.mincoverage,
                   snps=kwargs['snps'], pidiv_buffer=self.pidiv_buffer)
        
        return d

        
       
    
    
    
