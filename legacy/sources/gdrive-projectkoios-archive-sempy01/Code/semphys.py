from abc import ABC
from abc import abstractmethod
import numpy as np
import numpy.typing as npt

def _to_ndarray(
    x: npt.ArrayLike, 
    dtype=float
  ) -> np.ndarray:
  """
  Convert input to numpy array for consistent computation
  """
  return np.asarray(x, dtype=dtype)

class BandgapModel(ABC):
  """
  Abstract base class for temperature-dependent bandgap models.
  """

  @abstractmethod
  def Eg(self, T: npt.ArrayLike) -> np.ndarray:
    """Return the bandgap at temperature T [K]."""
    pass
  
class Varshi1970BandgapModel(BandgapModel):
  def __init__(self, Eg0: float, alpha: float, beta: float):
    self.Eg0 = Eg0
    self.alpha = alpha
    self.beta = beta

  def Eg(self, T: npt.ArrayLike) -> np.ndarray:
    T = _to_ndarray(T)
    return self.Eg0 - (self.alpha * T**2) / (T + self.beta)


class BoseEinsteinEquation(BandgapModel):
  def __init__(self, Eg0: float, aB: float, Theta: float):
    self.Eg0 = Eg0
    self.aB = aB
    self.Theta = Theta

  def Eg(self, T: npt.ArrayLike) -> np.ndarray:
    T = _to_ndarray(T)
    return self.Eg0 - self.aB / (np.exp(self.Theta / T) - 1)
