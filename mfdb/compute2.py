#################################################################
#                                                               #
# Database of Classical GL2 Holomorphic Modular Forms over Q    #
# (c) William Stein, 2012                                       #
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

# Basis Sage computations
@cached_function
def characters(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return [X[0] for X in DirichletGroup(N).galois_orbits()]

def character_to_int(eps):
    """
    Return integer corresponding to given character.
    """
    if eps.is_trivial():
        return 0
    N = eps.modulus()
    X = characters(N)
    try:
        return X.index(eps)
    except IndexError:
        # very unlikely -- would have to be some weird character
        # not got from character(N,i)
        for i, Y in enumerate(DirichletGroup(N).galois_orbits()):
            if X in Y:
                return i
        raise RuntimeError

#####################################################################
#
# Templates for the types of data we store. 
#
#####################################################################

class ObjTypeBase(object):
    def __init__(self, name, desc, params):
        self.name = name
        self.desc = desc
        self.params = params

    def __repr__(self):
        return 'Object type %s'%self.name
        
    def to_simple(self, obj):
        return self._to_simple(obj)

    def _to_simple(self, obj):
        return obj, 0

    def from_simple(self, obj, version):
        a = self._from_simple(obj, version)
        if a is None:
            raise ValueError, "version (=%s) not supported"%version
        return a
        
    def _from_simple(self, obj, version):
        if version == 0:
            return obj
        raise ValueError

#########################################################################
#
# The database supports storing the following object types.  Each
# class must derive from ObjTypeBase (as illustrated), and may provide
# a _to_simple and _from_simple method, which is used when saving the
# objects in the database or recovering them.  The _from_simple method
# should return None in case the given version is not supported.
#
#########################################################################

class ObjType_character(ObjTypeBase):
    def __init__(self, name):
        ObjTypeBase.__init__(self, name,
            'i-th representative of the Galois conjugacy classes of Dirichlet characters with modulus N',
            {'N':int,'i':int})
        
    def _to_simple(self, obj):
        version = 0
        return {'N':obj.modulus(), 'vals':list(obj.values_on_gens())}, version

    def _from_simple(self, obj, version):
        if version == 0:
            vals = obj['vals']
            return DirichletGroup(obj['N'], vals[0].parent())(vals)
        
            
class ObjType_ambient_modular_symbols_space(ObjTypeBase):
    def __init__(self, name):
        ObjTypeBase.__init__(self, name,
            'ambient modular symbols space of level N, weight k, and character i, with sign +1',
            {'N':int, 'k':int, 'i':int})

    def _to_simple(self, obj):
        version = 1
        i = character_to_int(obj.character())
        return {'space':(int(obj.level()), int(obj.weight()), int(i)),
                'eps':list(obj.character().values_on_gens()),
                'manin':[(t.i,t.u,t.v) for t in obj._manin_generators],
                'basis':obj._manin_basis,
                'rels':obj._manin_gens_to_basis,
                'mod2term':obj._mod2term}, version
    
    def _from_simple(self, obj, version):
        if version == 0:
            return obj
        elif version == 1:
            N     = obj['space'][0]
            k     = obj['space'][1]
            eps   = obj['eps']
            manin = obj['manin']
            basis = obj['basis']
            rels  = obj['rels']
            mod2term  = obj['mod2term']

            F = rels.base_ring()
            eps = DirichletGroup(N, F)(eps)
            
            from sage.modular.modsym.manin_symbols import ManinSymbolList, ManinSymbol
            manin_symbol_list = ManinSymbolList(k, manin)

            def custom_init(M):
                # reinitialize the list of Manin symbols with ours, which may be
                # ordered differently someday:
                syms = M.manin_symbols()
                ManinSymbolList.__init__(syms, k, manin_symbol_list)
                M._manin_generators = [ManinSymbol(syms, x) for x in manin]
                M._manin_basis = basis
                M._manin_gens_to_basis = rels
                M._mod2term = mod2term
                return M
            return ModularSymbols(eps, k, sign=1, custom_init=custom_init, use_cache=False)            

class ObjType_new_simple_modular_symbols_space(ObjTypeBase):
    def __init__(self, name):
        ObjTypeBase.__init__(self, name,
            'j-th new simple modular symbols symbols space of level N, weight k, and character i, with sign +1',
            {'N':int, 'k':int, 'i':int, 'j':int})

    def _to_simple(self, obj):
        version = 0
        obj.free_module().basis_matrix()
        
    def _from_simple(self, obj, version):


#######################

object_types = dict([(o.name, o) for o in
                     [o(k[len('ObjType_'):]) for k, o in
                                globals().iteritems() if k.startswith('ObjType_')]])

def ObjType(o):
    """
    Convert a string or object type o into an object type.
    """
    if isinstance(o, ObjTypeBase):
        return o
    return object_types[str(o)]

#####################################################################
#
# The database classes
#
#####################################################################


class Database(object):
    """
    Database of all computed modular forms information.

    All methods raise a KeyError if the requested data
    is not in the database.

    - hecke_eigenvalue_matrix(N, k, i, j, start, stop=None): matrix
      whose rows give the Hecke eigenvalues for primes p with start <=
      p < stop for the (N,k,i,j) new simple modular symbols space.

    - charpoly_list(N, k, i, j, start, stop=None): list whose entries
      are the characteristic polynomials of the Hecke eigenvalues a_p
      with start <= p < stop for the (N,k,i,j) new simple modular
      symbols space.

    - hecke_eigenvalue_field_basis(N, k, i, j): fixed choice of basis
      for the Hecke eigenvalue field for the (N,k,i,j) new simple
      modular symbols space.  The Hecke eigenvalues a_p are expressed
      in terms of this basis.

    - atkin_lehner_eigenvalues(N, k, i, j): sequence of integers +1/-1
      giving the eigenvalues of the Atkin-Lehner involutions w_{p^r}
      acting on the (N,k,i,j) new simple modular symbols space,
      ordered by p.
    
    - zeroes(N, k, i, j, max_imag=100): imaginary parts of zeros of
      L-series of conjugates of the newforms corresponding to the
      (N,k,i,j) new simple modular symbols space.  this is a list of
      lists of zeros of each conjugate.

    - lseries_leading(N, k, i, j): order of vanishing and leading
      coefficients of the L-series of conjugates of the newforms
      corresponding to the (N,k,i,j) new simple modular symbols space.
    """
    def __init__(self, to_simple=None, from_simple=None):
        if to_simple is None or from_simple is None:
            simple = SimpleMapper()
            to_simple = simple.obj_to_simple
            from_simple = simple.simple_to_obj
        self._obj_to_simple = to_simple
        self._simple_to_obj = simple.simple_to_obj

    def get(self, objtype, params):
        """
        Return the object of type objtype with given params stored in
        the database.  Raises a KeyError if there is no such object.
        """
        raise KeyError

    def set(self, objtype, params, value):
        """
        Store in the database the object value of type objtype defined
        by the given params.  
        """
        raise NotImplementedError

    def drop(self, objtype, confirm=False):
        """
        Delete everything in the database about the given objtype.
        """
        if not confirm:
            raise RuntimeError("are you sure?")
        objtype = ObjType(objtype)
        self._drop(objtype, confirm=confirm)

    def known(self, objtype):
        """
        Returns list of params of all objects with given objtype in
        the database.
        """
        raise NotImplementedError

class DatabaseFilesystem(Database):    
    def __init__(self, directory):
        self._directory = directory
        
class DatabaseNoSQLite(Database):
    RESERVED_COLUMNS = ['obj', 'version']
    
    def __init__(self, nsdb=None, **kwds):
        if nsdb is None:
            import nosqlite
            nsdb = nosqlite.Client('nosqlite_test_db').db
        self._nsdb = nsdb
        Database.__init__(self, **kwds)

    def _collection(self, objtype):
        objtype = ObjType(objtype)
        return self._nsdb.__getattr__(objtype.name)

    def get(self, objtype, **params):
        objtype = ObjType(objtype)        
        A = self._collection(objtype).find_one(**params)
        return objtype.from_simple(A['obj'], A['version'])
        
    def set(self, objtype, obj, **params):
        objtype = ObjType(objtype)        
        obj, version = objtype.to_simple(obj)
        C = self._collection(objtype)
        # changing before lets us store multiple versions
        # we could have a separate sweep to query and delete old versions
        params['version'] = version
        if C.count(**params) > 0:
            C.update({'obj':obj}, **params)
        else:
            # usual case
            params['obj'] = obj
            C.insert(**params)

    def _drop(self, objtype, confirm):
        if confirm:
            self._collection(objtype).delete()

    def known(self, objtype):
        objtype = ObjType(objtype) 
        C = self._collection(objtype)
        fields = [x for x in sorted(C._columns()) if x not in self.RESERVED_COLUMNS]
        return list(self._collection(objtype).find('', fields=fields))


class SimpleMapper(object):
    def obj_to_simple(self, objtype, obj):
        version = 0
        return obj, version
    
    def simple_to_obj(self, version, objtype, simple):
        return simple


################################################################
#                                                              #
# Structured mathematical interface to a database.             #
#                                                              #
################################################################

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
    

    
