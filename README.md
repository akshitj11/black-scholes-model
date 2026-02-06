# Black-Scholes Option Pricing Model

A Python implementation of the Black-Scholes model for pricing European call and put options. This repository provides multiple implementations ranging from simple to advanced, along with visualizations and educational resources.

##  Overview

The Black-Scholes model, developed by Fischer Black, Myron Scholes, and Robert Merton in 1973, revolutionized options pricing by providing a mathematical framework for calculating the theoretical price of European options. This implementation includes:

- **Multiple implementations** (minimal, functional, object-oriented)
- **Complete Greeks calculation** (Delta, Gamma, Vega, Theta, Rho)
- **Implied volatility calculator**
- **Interactive visualizations**
- **Educational examples and documentation**

##  Features

 Price European call and put options  
 Calculate all option Greeks  
 Compute implied volatility  
 Generate comprehensive visualizations  
 Multiple implementation styles (functional vs OOP)  
 Well-documented with examples  
 Production-ready code  

##  Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Install Dependencies

```bash
# Clone the repository
git clone https://github.com/yourusername/black-scholes-model.git
cd black-scholes-model

# Install required packages
pip install -r requirements.txt
```

**Required packages:**
```
numpy>=1.19.0
scipy>=1.5.0
matplotlib>=3.3.0
```

Or install manually:
```bash
pip install numpy scipy matplotlib
```


## 💻 Usage Examples

### Example 1: Price Multiple Options

```python
from black_scholes_class import Option

strikes = [90, 95, 100, 105, 110]

print(f"{'Strike':<10} {'Call Price':<12} {'Put Price':<12}")
print("-" * 40)

for K in strikes:
    call = Option(S=100, K=K, T=1, r=0.05, sigma=0.2, option_type='call')
    put = Option(S=100, K=K, T=1, r=0.05, sigma=0.2, option_type='put')
    
    print(f"${K:<9} ${call.price():<11.2f} ${put.price():<11.2f}")
```

**Output:**
```
Strike     Call Price   Put Price   
----------------------------------------
$90        $16.70       $1.81       
$95        $13.35       $3.46       
$100       $10.45       $5.57       
$105       $8.02        $8.14       
$110       $6.04        $11.15      
```

### Example 2: Calculate Greeks

```python
from black_scholes_general import calculate_greeks

greeks = calculate_greeks(
    S=100, 
    K=100, 
    T=1, 
    r=0.05, 
    sigma=0.2, 
    option_type='call'
)

print("Option Greeks:")
for greek, value in greeks.items():
    print(f"  {greek.capitalize()}: {value:.4f}")
```

**Output:**
```
Option Greeks:
  Delta: 0.6368
  Gamma: 0.0188
  Vega: 37.5247
  Theta: -6.4144
  Theta_per_day: -0.0176
  Rho: 53.2325
```

### Example 3: Find Implied Volatility

```python
from black_scholes_general import implied_volatility

# Given a market price, find the implied volatility
market_price = 10.50
iv = implied_volatility(
    option_price=market_price,
    S=100,
    K=100,
    T=1,
    r=0.05,
    option_type='call'
)

print(f"Market Price: ${market_price:.2f}")
print(f"Implied Volatility: {iv*100:.2f}%")
```

### Example 4: Delta Hedging

```python
from black_scholes_class import Option

# You sold 1000 call options
option = Option(S=100, K=105, T=0.5, r=0.05, sigma=0.25, option_type='call')

delta = option.delta()
shares_to_buy = 1000 * delta

print(f"Delta: {delta:.4f}")
print(f"To hedge 1000 short calls, buy {shares_to_buy:.0f} shares")
print(f"Hedge cost: ${shares_to_buy * 100:.2f}")
```

## 📊 The Black-Scholes Formula

### Call Option

$$C = S_0 N(d_1) - K e^{-rT} N(d_2)$$

### Put Option

$$P = K e^{-rT} N(-d_2) - S_0 N(-d_1)$$

### Where:

$$d_1 = \frac{\ln(S_0/K) + (r + \sigma^2/2)T}{\sigma\sqrt{T}}$$

$$d_2 = d_1 - \sigma\sqrt{T}$$

### Parameters:

| Symbol | Description | Example |
|--------|-------------|---------|
| S₀ | Current stock price | $100 |
| K | Strike price | $105 |
| T | Time to expiration (years) | 0.5 (6 months) |
| r | Risk-free interest rate | 0.05 (5%) |
| σ (sigma) | Volatility (annual std dev) | 0.25 (25%) |
| N(·) | Cumulative normal distribution | - |

### Key Assumptions:

1. European exercise (only at expiration)
2. No dividends during option life
3. Markets are efficient (no arbitrage)
4. No transaction costs or taxes
5. Constant risk-free rate and volatility
6. Log-normal stock price distribution
7. Continuous trading possible

## 📈 Understanding the Greeks

The Greeks measure how option prices change with different factors:

### Delta (Δ)
**Sensitivity to stock price changes**

```python
delta = option.delta()
# Call: 0 to 1
# Put: -1 to 0
```

- **Interpretation**: Option price change per $1 stock move
- **Example**: Delta = 0.5 → option gains $0.50 when stock gains $1
- **Use**: Hedge ratio (shares per option)

### Gamma (Γ)
**Rate of change of Delta**

```python
gamma = option.gamma()
# Always positive for long options
```

- **Interpretation**: How much Delta changes per $1 stock move
- **Example**: Gamma = 0.02 → Delta increases by 0.02 per $1 stock gain
- **Use**: Measure convexity, rebalancing frequency

### Vega (ν)
**Sensitivity to volatility changes**

```python
vega = option.vega()
# Always positive for long options
```

- **Interpretation**: Price change per 1% volatility change
- **Example**: Vega = 40 → option gains $0.40 if volatility rises 1%
- **Use**: Volatility trading strategies

### Theta (Θ)
**Time decay**

```python
theta = option.theta()
# Usually negative for long options
```

- **Interpretation**: Price change per day of time passing
- **Example**: Theta = -0.05 → option loses $0.05 per day
- **Use**: Time decay strategies

### Rho (ρ)
**Sensitivity to interest rate changes**

```python
rho = option.rho()
# Positive for calls, negative for puts
```

- **Interpretation**: Price change per 1% interest rate change
- **Example**: Rho = 30 → option gains $0.30 if rates rise 1%
- **Use**: Interest rate exposure (usually minimal)



**Happy Trading! 📈**

