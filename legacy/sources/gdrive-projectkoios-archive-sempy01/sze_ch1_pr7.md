![[Pasted image 20250930134108.png]]

**Temperature Dependence of the Bandgap**

The bandgap varation wr.t. temperature is described by the Varshni equation

##### Defintion: Varshni Equation
$$
E_g(T) = E_g(0) - \frac{\alpha T^2}{T + \beta}
$$
- $E_g$: bandgap, $\mathrm{eV}$
- $T$: temperature in $\mathrm{K}$
- $\beta$: fitting parameters
- $\alpha$: temperature coefficient, $\alpha$ (fitting parameters)

#### Parameters 
##### Silicon (Si)
- $E_g(0) = 1.17 \,\text{eV}$  
- $\alpha = 4.73 \times 10^{-4} \,\text{eV/K}$  
- $\beta = 636 \,\text{K}$  
##### Gallium Arsenide (GaAs)
- $E_g(0) = 1.519 \,\text{eV}$  
- $\alpha = 5.405 \times 10^{-4} \,\text{eV/K}$  
- $\beta = 204 \,\text{K}$  

##### Calculations

### Silicon
$$
E_g(100) = 1.17 - \frac{4.73 \times 10^{-4} \times 100^2}{100 + 636}
$$
$$
E_g(100) \approx 1.163 \,\text{eV}
$$
### GaAs
$$
E_g(100) = 1.519 - \frac{5.405 \times 10^{-4} \times 100^2}{100 + 204}
$$
$$
E_g(100) \approx 1.501 \,\text{eV}
$$

## Calculation at $T = 600 \,\text{K}$
### Silicon
$$
E_g(600) = 1.17 - \frac{4.73 \times 10^{-4} \times 600^2}{600 + 636}
$$
$$
E_g(600) \approx 1.036 \,\text{eV}
$$
### GaAs
$$
E_g(600) = 1.519 - \frac{5.405 \times 10^{-4} \times 600^2}{600 + 204}
$$
$$
E_g(600) \approx 1.350 \,\text{eV}
$$

## Summary

| Material | $E_g(100 \,\text{K})$       | $E_g(600 \,\text{K})$       |
| -------- | --------------------------- | --------------------------- |
| Si       | $\approx 1.163 \,\text{eV}$ | $\approx 1.036 \,\text{eV}$ |
| GaAs     | $\approx 1.501 \,\text{eV}$ | $\approx 1.350 \,\text{eV}$ |

##### Computationally
$$
E_g(T) = E_g(0) - \frac{\alpha T^2}{T + \beta}
$$
```python

# input: T: array
# parameters alpha, beta are parameters

def get_bandgap_varshni1970(
	T: float, 
	Eg_0: float,
	alpha: float, 
	beta: float) -> float:
  """
  T (float): temperature, in K
  Eg_0: bandgap at T=0
  alpha: fitting paramaeter
  beta: fitting parameter
  """
  
  E_g = Eg_0 - ((alpha*T**2)/(T+beta))
  return E_g


```

Varshni, Y. P. "Temperature dependence of the elastic constants." _Physical Review B_ 2.10 (1970): 3952.