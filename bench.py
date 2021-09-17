from .step import Step
from .ks import Ks

def get_bench(key):
  if key == 'step':
      bench = Step()
  elif key=='ks':
      bench = Ks()
  return bench
