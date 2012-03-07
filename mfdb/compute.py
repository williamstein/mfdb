import os

from sage.all import ModularSymbols, DirichletGroup, save, load, cputime, fork, parallel, Integer, version

class Filenames(object):
    def __init__(self, data):
        if not os.path.exists(data):
            raise RuntimeError, "please create the data directory '%s'"%data
        self._data = data
        
    def space(self, N, k, i):
        f = os.path.join(self._data, '%05d-%03d-%03d'%(N,k,i))
        if not os.path.exists(f):
            os.makedirs(f)
        return f

    def ambient(self, N, k, i):
        return os.path.join(self.space(N,k,i), 'M.sobj')

    def factor(self, N, k, i, d):
        f = os.path.join(self.space(N,k,i), '%03d'%d)
        if not os.path.exists(f):
            os.makedirs(f)
        return f
        
        
    def factor_basis_matrix(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'B.sobj')
    
    def factor_dual_basis_matrix(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'Bd.sobj')
    
    def factor_dual_eigenvector(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'v.sobj')

    def factor_eigen_nonzero(self, N, k, i, d):
        return os.path.join(self.factor(N,k,i,d), 'nz.sobj')

    def factor_aplist(self, N, k, i, d, *args):
        a = '-'.join('%05d'%x for x in args)
        return os.path.join(self.factor(N,k,i,d), 'aplist-%s.sobj'%a)

    def meta(self, filename):
        if filename.endswith('.sobj'):
            filename = filename[:-len('.sobj')]
        return filename + '-meta.sobj'

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
    
    eps = DirichletGroup(N).galois_orbits()[i][0]
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
    return load(filenames.ambient(N, k, i))


# decompositions

@fork    
def compute_decompositions(N, k, i):
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
    if not os.path.exists(filename):
        compute_ambient_space(N, k, i)
    if not os.path.exists(filename):
        return 
    
    eps = DirichletGroup(N).galois_orbits()[i][0]

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
    save(meta, filenames.meta(filename))

def compute_decomposition_ranges(Nrange, krange, irange, ncpu):
    @parallel(ncpu)
    def f(N,k,i):
        compute_decompositions(N,k,i)

    v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
    for X in f(v):
        print X
