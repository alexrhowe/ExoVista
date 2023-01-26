# distutils: language = c++
# distutils: sources = [Image.cpp]

cdef class PyImage:
    cdef Image *thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def SetupImage(self, rs, Te, rdb, ts, di, rd, drd, Qs):
        self.thisptr = new Image(rs, Te, rdb, ts, di, rd, drd, Qs)

    def disk_imager(self, x, y, z, r, dv, cosscattang):
        return self.thisptr.disk_imager(x, y, z, r, dv, cosscattang)