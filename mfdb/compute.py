import os

import sqlite3
import nosqlite

nsql = nosqlite.Client('nsql')

from sage.all import (ModularSymbols, DirichletGroup, trivial_character,
                      dimension_new_cusp_forms,
                      save, load,
                      cputime,
                      fork, parallel,
                      Integer,
                      prime_range, prime_divisors,
                      version,
                      Sequence,
                      cached_function,
                      Integer)

def dim_new(chi, k):
    if not isinstance(chi, (int, long, Integer)) and chi.is_trivial():
        return dimension_new_cusp_forms(chi.modulus(), k)
    else:
        return dimension_new_cusp_forms(chi, k)

@cached_function
def characters(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return [X[0] for X in DirichletGroup(N).galois_orbits()]

def character(N, i):
    if i==0:
        return trivial_character(N)
    print "character(%s, %s)"%(N,i)
    return characters(N)[i]

class Filenames(object):
    def __init__(self, data):
        if not os.path.exists(data):
            raise RuntimeError, "please create the data directory '%s'"%data
        self._data = data
        self._known_db_file = os.path.join(self._data, 'known.sqlite3')

    def space_name(self, N, k, i):
        return '%05d-%03d-%03d'%(N,k,i)
    
    def space(self, N, k, i, makedir=True):
        f = os.path.join(self._data, self.space_name(N,k,i))
        if makedir and not os.path.exists(f):
            os.makedirs(f)
        return f

    def ambient(self, N, k, i, makedir=True):
        return os.path.join(self.space(N,k,i,makedir=makedir), 'M.sobj')

    def decomp_meta(self, N, k, i):
        return os.path.join(self.space(N,k,i), 'decomp-meta.sobj')
    
    def factor(self, N, k, i, d, makedir=True):
        f = os.path.join(self.space(N,k,i,makedir), '%03d'%d)
        if makedir and not os.path.exists(f):
            os.makedirs(f)
        return f

    def number_of_known_factors(self, N, k, i):
        fn = self.space(N,k,i)
        return len([d for d in os.listdir(fn) if
                   d.isdigit() and os.path.isdir(os.path.join(fn, d))])
        
    def factor_basis_matrix(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'B.sobj')
    
    def factor_dual_basis_matrix(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'Bd.sobj')
    
    def factor_dual_eigenvector(self, N, k, i, d, makedir=True):
        return os.path.join(self.factor(N,k,i,d,makedir=makedir), 'v.sobj')

    def factor_eigen_nonzero(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'nz.sobj')

    def factor_aplist(self, N, k, i, d, makedir, *args):
        a = '-'.join('%05d'%x for x in args)
        return os.path.join(self.factor(N,k,i,d,makedir), 'aplist-%s.sobj'%a)

    def factor_atkin_lehner(self, N, k, i, d, makedir):
        return os.path.join(self.factor(N,k,i,d,makedir), 'atkin_lehner.txt')
    
    def meta(self, filename):
        if filename.endswith('.sobj'):
            filename = filename[:-len('.sobj')]
        return filename + '-meta.sobj'

    def find_known(self):
        """
        Return iterator of 5-tuples of Python ints, defined as follows:

            (N, k, i, newforms, maxp)
            (37, 2, 0, 2, 10000)

        Here N = level, k = weight, i = character, newforms = number of newforms,
        maxp = integer such that a_p is known for p<=maxp.

        If no newforms are known but there are newforms (they just
        haven't been computed), then newforms is set to -1.
        """
        for Nki in os.listdir(self._data):
            z = Nki.split('-')
            if len(z) == 3:
                N, k, i = [int(x) for x in z]
                newforms = [x for x in os.listdir(os.path.join(self._data, Nki)) if x.isdigit()]
                if len(newforms) == 0:
                    # maybe nothing computed?
                    if i == 0:
                        # program around a bug in dimension_new_cusp_forms: Trac 12640
                        d = dimension_new_cusp_forms(N)
                    else:
                        chi = character(N, i)
                        d = dimension_new_cusp_forms(chi, k)
                    if d == 0:
                        # definitely no newforms
                        yield (N,k,i,0,0)
                    else:
                        # we just don't know the newforms yet
                        yield (N,k,i,-1,0)
                else:
                    maxp = None
                    for n in newforms:
                        v = set([])
                        this_maxp = 0
                        for X in os.listdir(os.path.join(self._data, Nki, n)):
                            if X.startswith('aplist') and 'meta' not in X:
                                args = [int(a) for a in X.rstrip('.sobj').split('-')[1:]]
                                v.update(prime_range(*args))
                                this_maxp = max(this_maxp, max(args))
                        if len(v) != len(prime_range(this_maxp)):
                            # something missing!
                            print "data ranges are missing in the aplist data for %s"%Nki
                            maxp = 100
                        else:
                            maxp = this_maxp if maxp is None else min(this_maxp, maxp)

                    yield (N,k,i,len(newforms),maxp)

    def update_known_db(self):
        # 1. create the sqlite3 database
        if os.path.exists(self._known_db_file):
            os.unlink(self._known_db_file)
        db = sqlite3.connect(self._known_db_file)
        cursor = db.cursor()
        schema = 'CREATE TABLE "known" (N, k, i, newforms, maxp)'
        cursor.execute(schema)
        db.commit()

        # 2. iterate over known 5-tuples inserting them in db
        for t in self.find_known():
            print t
            cursor.execute("INSERT INTO known VALUES(?,?,?,?,?)", t)

        db.commit()


    def known(self, query):
        # TODO -- fix db
        query = query.replace('prec','maxp')
        
        # 1. open database
        db = sqlite3.connect(self._known_db_file)
        cursor = db.cursor()

        # 2. return result of query
        # really stupid and dangerous?
        if ';' in query:
            query = query.split(';')[0]

        cmd = 'SELECT * FROM known '
        if query.strip():
            cmd += 'WHERE %s'%query
        cmd += ' ORDER BY N, k, i'
        return cursor.execute(cmd)
        
    def find_missing(self, Nrange, krange, irange, fields=None):
        """
        Return generator of
        
             {'N':N, 'k':k, 'i':i, 'missing':missing, ...},

        where missing is in the intersection of fields and the
        following strings (or all strings if fields is None):

                'M', 'decomp',
                'aplist-00100',  'aplist-00100-01000',  'aplist-01000-10000',
                'charpoly-00100','charpoly-00100-01000','charpoly-01000-10000',
                'zeros', 
                'leading', 
                'atkin_lehner'

        If the string is not 'M' or 'decomp' then the meaning is that
        the indicated data is not complete for *all* newforms in this
        space, i.e., at least one is missing.  If 'decomp' is given,
        it means the decomp isn't complete (it could be partial).
        
        Spaces with dimension 0 are totally ignored. 

        Note that not every missing is listed, just the *next* one that
        needs to be computed, e.g., if (11,2,0,'M') is output, then
        nothing else for (11,2,0,*) will be output.
        """
        if fields is None:
            fields = set(['M', 'decomp',
                'aplist-00100',  'aplist-00100-01000',  'aplist-01000-10000',
                'charpoly-00100','charpoly-00100-01000','charpoly-01000-10000',
                'zeros', 
                'leading', 
                'atkin_lehner'])
        else:
            assert isinstance(fields, (list, tuple, str))
            fields = set(fields)

        space_params = set(os.listdir(self._data))
        for k in rangify(krange):
            for N in rangify(Nrange):
                for ch in rangify(irange):
                    
                    if isinstance(ch, str):
                        CHI = list(enumerate(characters(N)))
                        if ch == 'quadratic':
                            CHI = [(i,chi) for i,chi in CHI if chi.order()==2]
                        elif ch == 'all':
                            pass
                        else:
                            raise ValueError
                    else:
                        try:
                            CHI = [(ch, character(N, ch))]
                        except IndexError:
                            CHI = []
                            
                    for i, chi in CHI:
                        N,k,i = int(N), int(k), int(i)
                        obj = {'space':(N,k,i)}
                        dim = dim_new(chi, k)
                        if dim > 0:
                            fields0 = list(fields)
                            if 'M' in fields0:
                                fields0.remove('M')
                            if 'decomp' in fields0:
                                fields0.remove('decomp')
                            Nki = self.space_name(N,k,i)
                            d3 = os.path.join(self._data, Nki)
                            if Nki not in space_params or not os.path.exists(os.path.join(d3, 'M.sobj')):
                                if 'M' in fields:
                                    obj2 = dict(obj)
                                    obj2['missing'] = 'M'
                                    yield obj2
                                break
                            newforms = []
                            for fname in os.listdir(d3):
                                if fname.isdigit():
                                    # directory containing data about a newforms
                                    d2 = os.path.join(d3, fname)
                                    deg = os.path.join(d2, 'degree.txt')
                                    if os.path.exists(deg):
                                        degree = eval(open(deg).read())
                                    else:
                                        B_file = os.path.join(d2, 'B.sobj')
                                        if not os.path.exists(B_file):
                                            degree = None
                                        else:
                                            degree = load(B_file).nrows()
                                            open(os.path.join(d2, 'degree.txt'),'w').write(str(degree))
                                    f = {'fname':fname, 'degree':degree,
                                         'other':set([x.split('.')[0] for x in os.listdir(d2)])}
                                    newforms.append(f)
                            degs = [f['degree'] for f in newforms]
                            if None in degs:
                                sum_deg = None
                            else:
                                sum_deg = sum(degs)
                            if 'decomp' in fields and (sum_deg is None or sum_deg != dim):
                                obj2 = dict(obj)
                                obj2['missing'] = 'decomp'
                                obj2['newforms'] = newforms
                                if sum_deg > dim:
                                    obj2['bug'] = 'sum of degrees (=%s) is too big (should be %s) -- internal consistency error!'%(sum_deg, dim)
                                yield obj2
                                break
                            missing = []
                            for other in fields0:
                                for j in range(len(newforms)):
                                    if other not in newforms[j]['other']:
                                        missing.append(other)
                                        break
                            if missing:
                                missing.sort()
                                obj2 = dict(obj)
                                obj2['missing'] = 'other'
                                obj2['other'] = missing
                                yield obj2
                        
                            
################################################    

filenames = Filenames('data')

# ambient spaces

# use @fork to avoid any memory leaks
@fork    
def compute_ambient_space(N, k, i):
    if i == 'all':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn:
                compute_ambient_space(N,k,j)
        return

    if i == 'quadratic':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn and g[0].order()==2:
                compute_ambient_space(N,k,j)
        return

    filename = filenames.ambient(N, k, i)
    if os.path.exists(filename):
        return

    eps = character(N, i)
    t = cputime()
    M = ModularSymbols(eps, weight=k, sign=1)
    tm = cputime(t)
    save(M, filename)
    meta = {'cputime':tm, 'dim':M.dimension(), 'M':str(M), 'version':version()}
    save(meta, filenames.meta(filename))

def rangify(v):
    return [v] if isinstance(v, (int, long, Integer, str)) else v

def compute_ambient_spaces(Nrange, krange, irange, ncpu):
    @parallel(ncpu)
    def f(N,k,i):
        compute_ambient_space(N,k,i)

    v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
    for X in f(v):
        print X
    
    
def load_ambient_space(N, k, i):
    return load(filenames.ambient(N, k, i, makedir=False))

def load_factor(N, k, i, d, M=None):
    import sage.modular.modsym.subspace
    if M is None:
        M = load_ambient_space(N, k, i)
    f = filenames.factor(N, k, i, d, makedir=False)
    if not os.path.exists(f):
        raise RuntimeError, "no such factor (%s,%s,%s,%s)"%(N,k,i,d)
    B = load(filenames.factor_basis_matrix(N, k, i, d))
    Bd = load(filenames.factor_dual_basis_matrix(N, k, i, d))
    v = load(filenames.factor_dual_eigenvector(N, k, i, d))
    nz = load(filenames.factor_eigen_nonzero(N, k, i, d))
    B._cache['in_echelon_form'] = True
    Bd._cache['in_echelon_form'] = True
    # These two lines are scary, obviously, since they depend on
    # private attributes of modular symbols.
    A = sage.modular.modsym.subspace.ModularSymbolsSubspace(M, B.row_module(), Bd.row_module(), check=False)
    A._HeckeModule_free_module__dual_eigenvector = {('a',nz):(v,False)}
    A._is_simple = True
    A._HeckeModule_free_module__decomposition = {(None,True):Sequence([A], check=False)}
    return A


# decompositions

@fork    
def compute_decompositions(N, k, i):
    if i == 'all':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn:
                compute_decompositions(N,k,j)
        return

    if i == 'quadratic':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn and g[0].order()==2:
                compute_decompositions(N,k,j)
        return

    filename = filenames.ambient(N, k, i)
    if not os.path.exists(filename):
        print "Ambient space (%s,%s,%s) not computed."%(N,k,i)
        return
        #compute_ambient_space(N, k, i)
    if not os.path.exists(filename):
        return 
    
    eps = DirichletGroup(N).galois_orbits()[i][0]

    if os.path.exists(filenames.factor(N, k, i, 0, makedir=False)):
        return
    
    t = cputime()
    M = load_ambient_space(N, k, i)
    D = M.cuspidal_subspace().new_subspace().decomposition()
    for d in range(len(D)):
        f = filenames.factor_basis_matrix(N, k, i, d)
        if os.path.exists(f):
            continue
        A = D[d]
        B  = A.free_module().basis_matrix()
        Bd = A.dual_free_module().basis_matrix()
        v  = A.dual_eigenvector(names='a', lift=False)    # vector over number field
        nz = A._eigen_nonzero()
        
        save(B, filenames.factor_basis_matrix(N, k, i, d))
        save(Bd, filenames.factor_dual_basis_matrix(N, k, i, d))
        save(v, filenames.factor_dual_eigenvector(N, k, i, d))
        save(nz, filenames.factor_eigen_nonzero(N, k, i, d))
        
    tm = cputime(t)
    meta = {'cputime':tm, 'number':len(D), 'version':version()}
    save(meta, filenames.decomp_meta(N, k, i))

def compute_decomposition_ranges(Nrange, krange, irange, ncpu):
    @parallel(ncpu)
    def f(N,k,i):
        compute_decompositions(N,k,i)

    v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
    for X in f(v):
        print X

def atkin_lehner_signs(A):
    N = A.level()
    return [A._compute_atkin_lehner_matrix(p)[0,0] for p in prime_divisors(N)]

# atkin_lehner
@fork    
def compute_atkin_lehner(N, k, i):
    filename = filenames.ambient(N, k, i)
    if not os.path.exists(filename):
        print "Ambient (%s,%s,%s) space not computed."%(N,k,i)
        return
        #compute_ambient_space(N, k, i)

    print "computing atkin-lehner for (%s,%s,%s)"%(N,k,i)
    m = filenames.number_of_known_factors(N, k, i)    
    M = load_ambient_space(N, k, i)
    for d in range(m):
        atkin_lehner_file = filenames.factor_atkin_lehner(N, k, i, d, False)
        if os.path.exists(atkin_lehner_file):
            print "skipping computing atkin_lehner for (%s,%s,%s,%s) since it already exists"%(N,k,i,d)
            # already done
            continue
        
        # compute atkin_lehner
        print "computing atkin_lehner for (%s,%s,%s,%s)"%(N,k,i,d)
        t = cputime()
        A = load_factor(N, k, i, d, M)
        al = ' '.join(['+' if a > 0 else '-' for a in atkin_lehner_signs(A)])
        print al
        open(atkin_lehner_file, 'w').write(al)
        tm = cputime(t)
        meta = {'cputime':tm, 'version':version()}
        save(meta, filenames.meta('atkin_lehner'))


# aplists
    
@fork    
def compute_aplists(N, k, i, *args):
    if i == 'all':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn:
                compute_aplists(N,k,j,*args)
        return

    if i == 'quadratic':
        G = DirichletGroup(N).galois_orbits()
        sgn = (-1)**k
        for j, g in enumerate(G):
            if g[0](-1) == sgn and g[0].order()==2:
                compute_aplists(N,k,j,*args)
        return

    if len(args) == 0:
        args = (100, )

    filename = filenames.ambient(N, k, i)
    if not os.path.exists(filename):
        print "Ambient (%s,%s,%s) space not computed."%(N,k,i)
        return
        #compute_ambient_space(N, k, i)

    print "computing aplists for (%s,%s,%s)"%(N,k,i)
        
    m = filenames.number_of_known_factors(N, k, i)

    if m == 0:
        # nothing to do
        return
    
    M = load_ambient_space(N, k, i)
    for d in range(m):
        aplist_file = filenames.factor_aplist(N, k, i, d, False, *args)
        if os.path.exists(aplist_file):
            print "skipping computing aplist(%s) for (%s,%s,%s,%s) since it already exists"%(args, N,k,i,d)
            # already done
            continue
        
        # compute aplist
        print "computing aplist(%s) for (%s,%s,%s,%s)"%(args, N,k,i,d)
        t = cputime()
        A = load_factor(N, k, i, d, M)
        aplist, _ = A.compact_system_of_eigenvalues(prime_range(*args), 'a')
        print aplist, aplist_file
        save(aplist, aplist_file)
        tm = cputime(t)
        meta = {'cputime':tm, 'version':version()}
        save(meta, filenames.meta(aplist_file))

def compute_aplists_ranges(Nrange, krange, irange, ncpu, *args):
    @parallel(ncpu)
    def f(N,k,i):
        compute_aplists(N,k,i,*args)

    v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
    for X in f(v):
        print X

def phase1_goals(stages, fields=None):
    """
      Trivial character S_k(G0(N)):
         0. k=2   and N<=4096= 2^12
         1. k<=16 and N<=512 = 2^9
         2. k<=32 and N<=32  = 2^5
      Nontrivial character S_k(N, chi):
         3. k<=16, N<=128 = 2^7, all chi quadratic
         4. k=2,   N<=128 = 2^7, all chi!=1, up to Galois orbit
    """
    if 0 in stages:
        for X in filenames.find_missing(range(1,4096+1), [2], 0, fields=fields):
            yield X
            
    if 1 in stages:
        for X in filenames.find_missing(range(1,513),
                                        range(2,17,2), 0, fields=fields):
            yield X
            
    if 2 in stages:
        for X in filenames.find_missing(range(1,33),
                                        range(2,33,2), 0, fields=fields):
            yield X
            
    if 3 in stages:
        for X in filenames.find_missing(range(1,129),
                                 range(2,17), 'quadratic', fields=fields):
            yield X
            
    if 4 in stages:
        for X in filenames.find_missing(range(1,129), 2, 'all', fields=fields):
            yield X
        

def phase1_goals0(stages, fields=None):
    v = []
    for X in phase1_goals(stages, fields):
        print X
        v.append(X)
    return v

# suggestion for collection is: mfdb.compute.nsql.missing.phase0
def phase1_goals_db(collection, stages, fields=None):
    for X in phase1_goals(stages, fields):
        print X
        collection.insert(X)

def generate_computations_missing_M(collection):
    v = []
    for X in collection.find(missing="M"):
        N,k,i = X['space']
        v.append('compute_ambient_space(%s,%s,%s)'%(N,k,i))
    return list(sorted(set(v)))
        
def generate_computations_missing_decomp(collection):
    v = []
    for X in collection.find(missing="decomp"):
        N,k,i = X['space']
        v.append('compute_decompositions(%s,%s,%s)'%(N,k,i))
    return list(sorted(set(v)))        
        
def generate_computations_missing_aplists(collection):
    v = []
    for X in collection.find(missing="other"):
        N,k,i = X['space']
        if 'aplist-00100' in X['other']:
            v.append('compute_aplists(%s, %s, %s, 100)'%(N,k,i))
        if 'aplist-00100-01000' in X['other']:
            v.append('compute_aplists(%s, %s, %s, 100, 1000)'%(N,k,i))
        if 'aplist-01000-10000' in X['other']:
            v.append('compute_aplists(%s, %s, %s, 1000, 10000)'%(N,k,i))
##         if 'charpoly-00100' in X['other']:
##             v.append('compute_charpolys(%s, %s, %s, 100)'%(N,k,i))
##         if 'charpoly-00100-01000' in X['other']:
##             v.append('compute_charpolys(%s, %s, %s, 100, 1000)'%(N,k,i))
##         if 'charpoly-01000-10000' in X['other']:
##             v.append('compute_charpolys(%s, %s, %s, 1000, 10000)'%(N,k,i))
    return list(sorted(set(v)))        
    


import sage.all
def parallel_eval(v, ncpu, do_fork=True):
    """
    INPUT:
    - v -- list of strings that can be eval'd in this scope.
    """
    @parallel(ncpu)
    def f(X):
        if do_fork:
            @fork
            def g():
                return eval(X)
            return g()
        else:
            return eval(X)
    for Z in f(v):
        yield Z

        
        
        
    
