from sage.all import (DirichletGroup,
                      )

class Service(object):
    def __init__(self, db_file):
        self._db_file = db_file

    def update_database(self, filenames):
        """
        Create the sqlite3 database using the data pointed to by the
        filenames object.  This database has schema:

            N   k  i   newforms   maxp
            37  2  0   2          2

            If no newforms have been computed, that column of the
            table is empty.  The above can all be deduced from the
            filesystem, without having to open any sobj files.
        """
        raise NotImplementedError
        
    def characters(self, N):
        """
        Return a sorted list of characters divided by Galois orbit,
        with their sequence number and parity.
        """
        return DirichletGroup(N).galois_orbits()
    
    def known(self, query):
        """
        Given a sqlite query string involving N,k,i,newforms,maxp,
        return all corresponding 5-tuples:
        
            (level, weight, character, number of newforms,
                   largest n such that a_p is known for p<=maxp)

        If no newforms have yet been computed, but the space is known,
        then the number of newforms is set to -1.
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
