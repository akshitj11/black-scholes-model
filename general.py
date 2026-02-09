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
     Calculate all Greeks for an option
    
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
    option_type : str
        'call' or 'put'
    
    Returns:
    --------
    dict : Dictionary with all Greeks
    
    Example:
    --------
    >>> greeks = calculate_greeks(100, 100, 1, 0.05, 0.2, 'call')
    >>> print(f"Delta: {greeks['delta']:.4f}")
    """
    # Calculate d1 and d2
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    #delta
    if option_type == 'call':
        delta = norm.edf(d1)
    else:
        delta = norm.edf(d1) - 1

    #gamma 
    gamma = norm.pdf(d1)/(S*sigma*np.sqrt(T))

    #vega
    vega = S*norm.pdf(d1) * np.sqrt(T)

    #theta
    first_term = -(S*norm.pdf(d1)*sigma) / (2*np.sqrt(T))
    if option_type == 'call':
        second_term= -r*K*np.exp(-r*T)*norm.cdf(d2)
    else:
        second_term=r*K*np.exp(-r*T) * norm.cdf(-d2)
    theta = first_term + second_term

    #rho
    if option_type == 'call':
        rho = K*T*np.exp(-r*T)*norm.cdf(d2)
    else:
        rho = -K *T*np.exp(-r*T)*norm.cdf(d2)

    return{
        'delta': delta ,
        'gamma': gamma ,
        'vega': vega,
        'theta':theta ,
        'theta_per_day': theta / 365 ,
        'rho': rho
    }    




    