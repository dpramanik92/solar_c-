from ctypes import cdll
import ctypes
import numpy as np



lib = cdll.LoadLibrary('./libchi2.so')

class solchi2(object):
    def __init__(self):
        print(" Initializing ... ")
        
    def Set_resolution(self,what,samp_min,samp_max,samp_num):
        lib.solSetRes.argtypes = [ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_int]
        if(what=='SOL_YES'):
            lib.solSetRes(1,samp_min,samp_max,samp_num)
            print("Resolution function set...")
        if(what=='SOL_NO'):
            lib.solSetRes(0,samp_min,samp_max,samp_num)
            print("Resolution function wiil not be taken")
        
    
    def init_exp(self,experiment,which_type,channel):

        lib.solInitExp.argtypes = [ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int]
        lib.solInitExp(str(experiment).encode('ascii'),str(which_type).encode('ascii'),channel)

    def calChi2sing(self,params):
        lib.solChi2.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        lib.solChi2.restype = ctypes.c_double
        f = lib.solChi2(params)
        return f


    def calChi2(self,params):
        X = []
        for i in range(0,len(params)):
            X.append(params[i])

        N = len(X[0])
        X = np.array(X)

        lib.solChi2.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        lib.solChi2.restype = ctypes.c_double

        r = []
        for i in range(0,N):
            A = X.T[i]
            z = lib.solChi2(np.array(tuple(A)))
            r.append(z)

        r = np.array(r)

        return r




