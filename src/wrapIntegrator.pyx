# distutils: language = c++
# distutils: sources = [Integrator.cpp]

cdef class PyIntegrator:
    cdef Integrator *thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def SetupIntegrator(self, GMin, stin, tin, dtin):
        self.thisptr = new Integrator(GMin, stin, tin, dtin)

    def integrate(self):
        return self.thisptr.integrate()

    def getState(self):
        return self.thisptr.getState()

    def getTime(self):
        return self.thisptr.getTime()

    def getDeltaT(self):
        return self.thisptr.getDeltaT()

    def setDeltaT(self, newdeltat):
        return self.thisptr.setDeltaT(newdeltat)

    def equations_of_motion(self, local_state):
        return self.thisptr.equations_of_motion(local_state)