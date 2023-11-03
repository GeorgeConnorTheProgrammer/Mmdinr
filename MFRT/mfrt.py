import ctypes
import sys
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))

libMFRT = ctypes.cdll.LoadLibrary(dir_path + "/mfrt.dll")

libMFRT.run_model.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int]
libMFRT.run_model(0.1, 50, 300, 11, 12)