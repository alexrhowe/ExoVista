from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Integrator.h":
    cdef cppclass Integrator:
        Integrator(vector[double] GMin, vector[double] stin, double tin, double dtin) except +
        void integrate()
        vector[double] getState()
        double getTime()
        double getDeltaT()
        void setDeltaT(double newdeltat)
        vector[double] equations_of_motion(vector[double] local_state)