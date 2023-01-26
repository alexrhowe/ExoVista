from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Image.h":
    cdef cppclass Image:
        Image(double rs, double Te, double rdb, double ts, vector[vector[double]] di, vector[double] rd, vector[double] drd, vector[vector[double]] Qs) except +
        vector[vector[vector[double]]] disk_imager(vector[vector[vector[double]]] x, vector[vector[vector[double]]] y, vector[vector[vector[double]]] z, vector[vector[vector[double]]] r, vector[vector[vector[double]]] dv, vector[vector[vector[double]]] cosscattang)