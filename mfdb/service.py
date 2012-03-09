import compute

class Service(object):
    def __init__(self, filenames):
        self._filenames = filenames

    def update_database(self):
        """
        Create the sqlite3 database using the data pointed to by the
        filenames object.  This database has schema:

            N   k  i   newforms   maxp
            37  2  0   2          2

            If no newforms have been computed, that column of the
            table is empty.  The above can all be deduced from the
            filesystem, without having to open any sobj files.
        """
        self._filenames.update_known_db()
        
    def characters(self, N):
        """
        Return a sorted list of characters divided by Galois orbit,
        with their sequence number and parity.
        """
        return compute.characters(N)
    
    def known(self, query):
        """
        Given a sqlite query string involving N,k,i,newforms,maxp,
        return all corresponding 5-tuples:
        
            (level, weight, character, number of newforms,
                   largest n such that a_p is known for p<=maxp)

        If no newforms have yet been computed, but the space is known,
        then the number of newforms is set to -1.
        """
        return self._filenames.known(query)

    def modsym(self, N, k, i, j=None):
        """
        Return the corresponding simple modular symbols factor.
        """
        if j is None:
            return compute.load_ambient_space(N, k, i)
        else:
            return compute.load_factor(N, k, i, j)

    def aplist(self, N, k, i, j, pmax):
        """
        Return the pair (v, aplist), where v is the dual eigenvector,
        aplist is a matrix and v*aplist is the list of coefficients in
        terms of a power basis.
        """
        raise NotImplementedError
        #v = load(self._filenames.factor_dual_eigenvector(N, k, i, d))
        #    files = self._filenames.factor_aplist_pmax(N, k, i, d, pmax)
