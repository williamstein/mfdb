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
    def __init__(self):
        raise NotImplementedError, "derive!"

    def ambient_modular_symbols_space(self, N, k, i):
        raise NotImplementedError

    def newform_modular_symbols_space(self, N, k, i, j):
        raise NotImplementedError

    def newform_aplist(self, N, k, i, j, start, stop=None):
        raise NotImplementedError

    def newform_atkin_lehner_eigenvalues(self, N, k, i, j):
        raise NotImplementedError
    
class DatabaseFilesystem(Database):    
    def __init__(self, directory):
        self._directory = directory
        
class DatabaseNoSQLite(Database):    
    def __init__(self, nsdb):
        self._nsdb = nsdb

class AmbientSpaces(object):
    """
    All ambient spaces in the database.

        sage: X = AmbientSpaces()
        sage: M = X[389,2,0]
        sage: F = M.newform_classes()
        sage: c = F[0]
        sage: c.atkin_lehners()
        sage: c.degree()
        sage: c.hecke_eigenvalue_field()
        sage: c.q_expansions(10)
        sage: f = c.newforms()[0]
        sage: f.q_expansion(10)
        sage: L = f.lseries()
        sage: L.zeros()
        sage: L.lseries_leading()
    """
    def __call__(self, N, k, i):
        return AmbientSpace(N, k, i)
    
    def missing(self, Nrange, krange, irange):
        """
        Return list of parameters (N,k,i) such that N in Nrange, k in
        krange, and i in irange, such that the space (N,k,i) has not
        been computed.  irange can be 0, 'all', 'quadratic', or a list
        of integers.
        """
        raise NotImplementedError

    def known(self):
        """
        Return list of parameters (N,k,i) such that a presentation for
        the ambient space of modular symbols of level N, weight k, and
        character defined by i is in the database.
        """
        raise NotImplementedError


class NewformClasses(object):
    """
    All newform classes in the database.
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
        """
        Return the newform with parameters (N,k,i,j), where N is a
        positive integer, k an integer >= 2, and i and j are
        nonnegative integers.
        """
        return AmbientSpace(N,k,i)[j]

    def properties(self):
        """
        Return list of the properties whose computation is supported.
        """
        list(sorted(self._properties.keys()))

    def missing(self, Nrange, krange, irange, properties=None):
        """
        Return dictionary {(N,k,i,j):['prop1', 'prop2', ...], ...} where
          N in Nrange, k in krange, and i in irange
        and the 'propi' are the unknown properties for some newform in
        the space.  Here properties is a sublist of self.properties()
        """
        raise NotImplementedError

    def known(self, properties=None):
        """
        Return list of 4-tuples (N,k,i,j) of all newforms in the
        database that have all given properties computed.
        """
        raise NotImplementedError

class AmbientSpace(object):
    """
    An ambient space of modular forms.
    """
    def __init__(self, N, k, i):
        self._params = (N, k, i)        

    def __repr__(self):
        return "Ambient space (%s,%s,%s)"%self._params

    def __len__(self):
        """
        Return number of Galois conjugacy classes of newforms in this space.
        """
        raise NotImplementedError

    def __getitem__(self, j):
        """
        Return the j-th Galois conjugacy class of newforms in this space.
        """
        if j < 0:
            j += len(self)
        if j < 0 or j >= len(self): raise IndexError
        N,k,i = self._params
        return Newform(N,k,i,j)

    def newform_classes(self):
        """
        Return all Galois conjugacy classes of newforms in this space.
        """
        return NewformClasses()(*self._params)

class NewformClass(object):
    """
    A Gal(Qbar/Q) orbit of newforms.
    """
    def __init__(self, N, k, i, j):
        self._params = (N,k,i,j)

    def __repr__(self):
        return "Newform (%s,%s,%s,%s)"%self._params

    def newforms(self):
        """
        Return all the newforms in this class.
        """
        raise NotImplementedError

    def __getitem__(self, i):
        raise NotImplementedError
    
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
        """
        Compute and store zeros up to the given bound for all the
        conjugate newforms in this class.
        """
        raise NotImplementedError
    
    ###################################################################

    def degree(self):
        """
        Return the degree of this newform, which is the degree of the
        Hecke eigenvalue field.
        """
        raise NotImplementedError

    def atkin_lehner_eigenvalues(self):
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
        the algebraic integers a_p, indexed by primes p<B.
        """
        raise NotImplementedError

    def hecke_eigenvalue_field(self):
        """
        Return abstract number field generated by the coefficients of
        this newform.
        """
        raise NotImplementedError

    def q_expansion(self, B=infinity):
        """
        Return q-expansion up `O(q^B)` with coefficients in the
        abstract Hecke eigenvalue field, computed using as many
        coefficients as we know.
        """
        raise NotImplementedError

class Newform(object):
    """
    A newform equipped with a specific choice of embedding of its
    coefficients into the complex numbers.
    """
    def __init__(self, newform_class, n):
        self._newform_class = newform_class
        self._n = n
        
    def hecke_eigenvalues(self, B=infinity, prec=53):
        """
        Return dictionary of all known Hecke eigenvalues a_p, indexed
        by primes p<B.
        """
        raise NotImplementedError

    def q_expansion(self, B=infinity, prec=53):
        """
        Return q-expansion up to `O(q^B)` with coefficients in complex
        numbers, computed using as many coefficients as we know.
        """
        raise NotImplementedError

    def lseries(self):
        """
        Return the L-series of this newform.
        """
        return LSeries(self)

class LSeries(object):
    """
    The L-series attached to a newform.
    """
    def __init__(self, newform):
        self._newform = newform
    
    def zeros(self, B=infinity):
        """
        Return imaginary parts of known zeros of the L-function of
        this newform on the center of the critical strip with
        imaginary part bounded by B in absolute value.
        """
        raise NotImplementedError

    def lseries_leading(self):
        """
        Return dictionary {i:(ord, leading), ...} with keys all the
        integers inside the critical strip and values
        
            ord     = order of vanishing, and
            leading = leading coefficient
            
        of this L-function.
        """
        raise NotImplementedError
    

    
