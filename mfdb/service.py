class Service(object):
    def characters(self, N):
        """
        Return a sorted list of characters divided by Galois orbit,
        with their sequence number and parity.
        """
        raise NotImplementedError
    
    def known_newforms(self, query):
        """
        Given a sqlite query string involving N,k,i, return all spaces
        (N,k,i) for which a decomposition has been computed, where
        N,k,i satisfy the given constraint.
        """
        raise NotImplementedError

    def known_aplists(self, query):
        """
        Given a sqlite query string involving N,k,i,maxp return all
        spaces (N,k,i) for which a decomposition has been computed,
        where N,k,i,maxp satisfy the given constraint.
        """
        raise NotImplementedError

    def count(self, query):
        """
        Return the number of newforms in the spaces (N,k,i) satisfying
        the given query string involving N,k,i.
        """
        raise NotImplementedError

    def modsym(self, N, k, i, j):
        """
        Return the corresponding simple modular symbols factor.
        """
        raise NotImplementedError

    def aplist(self, N, k, i, j, pmax):
        """
        Return the pair (v, aplist), where v is the dual eigenvector,
        aplist is a matrix and v*aplist is the list of coefficients in
        terms of a power basis.
        """
        raise NotImplementedError
