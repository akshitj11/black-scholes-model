"""
implementation for pricing european options 
"""

import numpy as np
from scipy.stats import norm 

def black_scholes_call(S,K,T,r,sigma):
    """
    Calculate Black-Scholes price for a European call option
    
    Parameters:
    -----------
    S : float
        Current stock price
    K : float
        Strike price
    T : float
        Time to expiration (in years)
    r : float
        Risk-free interest rate (annual, as decimal)
    sigma : float
        Volatility (annual, as decimal)
    
    Returns:
    --------
    float : Call option price
    
    Example:
    --------
    >>> price = black_scholes_call(100, 100, 1, 0.05, 0.2)
    >>> print(f"Call price: ${price:.2f}")
    """
    # Calculate d1 and d2
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    #calculate call price
    call_price = S * norm.cdf(d1) - K * np.exp(-r*T) * norm.cdf(d2)

    return call_price 

def black_scholes_put(S,K,T,r,sigma):
    """
    Calculate Black-Scholes price for a European put option
    
    Parameters:
    -----------
    S : float
        Current stock price
    K : float
        Strike price
    T : float
        Time to expiration (in years)
    r : float
        Risk-free interest rate (annual, as decimal)
    sigma : float
        Volatility (annual, as decimal)
    
    Returns:
    --------
    float : Put option price
    
    Example:
    --------
    >>> price = black_scholes_put(100, 100, 1, 0.05, 0.2)
    >>> print(f"Put price: ${price:.2f}")
    """

     # Calculate d1 and d2
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    #calculate put price 
    put_price = K * np.exp(-r*T) * norm.cdf(-d2) - S * norm.cdf(-d1)
  
    return put_price 

def calculate_greeks(S,K,T,r,sigma ,option_type='call'):
    """
    """
    