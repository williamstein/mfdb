#################################################################
#                                                               #
# Database of Classical GL2 Holomorphic Modular Forms over Q    #
#                                                               #
#################################################################

# system imports
import os

# database imports
import nosqlite

# sage imports
from sage.all import (ModularSymbols,
                      DirichletGroup, trivial_character,
                      dimension_new_cusp_forms,
                      save, load,
                      version,
                      cputime,
                      fork, parallel, cached_function,
                      Integer, infinity,
                      prime_range, prime_divisors,
                      Sequence)


# The database classes

class Database(object):
    """
    Database of all computed modular forms information.
    """
    def __init__(self, db):
        self._db = db
        
class AmbientSpaces(object):
    """
    All ambient spaces in the database.
    """
    def __call__(self, N, k, i):
        return AmbientSpace(N,k,i)
    
    def missing(self, Nrange, krange, irange):
        raise NotImplementedError

    def known(self):
        raise NotImplementedError

class Newforms(object):
    """
    All newforms in the database.
    """
    _properties = {
        'atkin_lehner':{'depends':[], 'command':('compute_atkin_lehner',)},

        'aplist-00100':{'depends':[], 'command':('compute_aplist', 100)},
        'aplist-00100-01000':{'depends':[], 'command':('compute_aplist', 100, 1000)},
        'aplist-01000-10000':{'depends':[], 'command':('compute_aplist', 1000, 10000)},

        'charpoly-00100':{'depends':['aplist-00100'], 'command':('compute_charpoly', 100)},
        'charpoly-00100-01000':{'depends':['aplist-00100-01000'], 'command':('compute_charpoly', 100, 1000)},
        'charpoly-01000-10000':{'depends':['aplist-01000-10000'], 'command':('compute_charpoly', 1000, 10000)},

        'leading':{'depends':['aplist-00100','aplist-00100-01000','aplist-01000-10000'],
                       'command':('compute_leading')},

        'zeros-00100':{'depends':['aplist-00100','aplist-00100-01000','aplist-01000-10000'],
                       'command':('compute_zeros', 100)},
    }
    
    def __call__(self, N, k, i, j):
        return AmbientSpace(N,k,i)[j]

    def properties(self):
        list(sorted(self._properties.keys()))

    def missing(self, Nrange, krange, irange, properties=None):
        """
        Return dictionary {(N,k,i,j):['prop1', 'prop2', ...], ...} where
          N in Nrange, k in krange, and i in irange
        and the 'propi' are the unknown properties for some newform in
        the space.  Here properties is a sublist of self.properties()
        """
        raise NotImplementedError

    def known(self):
        raise NotImplementedError

class AmbientSpace(object):
    """
    A specific ambient space of modular forms.
    """
    def __init__(self, N, k, i):
        self._params = (N,k,i)        

    def __repr__(self):
        return "Ambient space (%s,%s,%s,%s)"%self._params

    def __len__(self):
        raise NotImplementedError

    def __getitem__(self, j):
        if j < 0:
            j += len(self)
        if j < 0 or j >= len(self): raise IndexError
        N,k,i = self._params
        return Newform(N,k,i,j)

    def newforms(self):
        return Newforms()(*self._params)

class Newform(object):
    """
    A specific newform.
    """
    def __init__(self, N, k, i, j):
        self._params = (N,k,i,j)

    def __repr__(self):
        return "Newform (%s,%s,%s,%s)"%self._params

    ###################################################################
    
    def compute_atkin_lehner(self):
        """
        Compute and store in the database the Atkin-Lehner signs.
        """
        raise NotImplementedError

    def compute_hecke_eigenvalues(self, start, stop=None):
        """
        Compute and store in the database the Hecke eigenvalues a_p
        for p in prime_range(start, stop).
        """
        raise NotImplementedError

    def compute_charpoly(self, start, stop=None):
        """
        Compute and store in the database the characteristic polynomial
        of the Hecke eigenvalue a_p for p in prime_range(start, stop).
        """
        raise NotImplementedError

    def compute_leading(self):
        """
        Compute and store leading coefficient of the L-series.
        """
        raise NotImplementedError

    def compute_zeros(self, max_imag):
        raise NotImplementedError


    ###################################################################

    def atkin_lehners(self):
        """
        Return list of integers +1 or -1 that are the eigenvalues of
        the Atkin-Lehner involutions acting on this newform,
        corresponding to the prime-power divisors of the level (in
        order of prime).
        """
        raise NotImplementedError

    def charpolys(self, B=infinity):
        """
        Return dictionary of all known characteristic polynomials of
        the a_p, indexed by primes p<B.
        """
        raise NotImplementedError

    def hecke_eigenvalues(self, embedding=None, B=infinity):
        """
        Return dictionary of all known Hecke eigenvalues a_p, indexed
        by primes p<B.
        """
        raise NotImplementedError

    def zeros(self, B=infinity):
        """
        Return imaginary parts of known zeros of the L-function on the
        center of the critical strip with imaginary part bounded by B
        in absolute value.
        """
        raise NotImplementedError

    def lseries_leading(self):
        """
        Return dictionary {i:(ord, leading), ...} with keys the
        integers inside the critical strip and values
            ord     = order of vanishing, and
            leading = leading coefficient.
        """
        raise NotImplementedError
    
    ###################################################################

    def q_expansion(self, embedding=None, B=infinity):
        """
        Return q-expansion computed using as many coefficients as we
        know.
        """
        raise NotImplementedError

    def lseries(self):
        """
        Return L-series computed using as many coefficients as we know.
        """
        raise NotImplementedError

    

    
