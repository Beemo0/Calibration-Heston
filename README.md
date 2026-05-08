# Heston Model Calibration and Volatility Smile
 
## Overview
 
This project implements the **Heston stochastic volatility model** for European option pricing via **Monte-Carlo simulation**, calibrates the model parameters using the **Levenberg-Marquardt algorithm**, and reconstructs the **implied volatility smile** via the **Newton algorithm**.
 
Unlike Black-Scholes which assumes constant volatility, the Heston model treats volatility as a stochastic process with mean-reversion. This allows it to capture market phenomena such as **volatility clustering**, **skewness**, and **fat tails** in the distribution of returns.
 
The implementation is written in **Python** and relies on Monte-Carlo simulation with variance reduction techniques.
 
---
 
## Project Structure
 
```
TP5/
├── TP5.py                    # Python implementation of all computations
└── README.md                 # Project documentation
```
 
---
 
# Theoretical Background
 
## Heston Model Dynamics
 
The asset price $S_t$ and its variance $v_t$ follow the coupled SDEs:
 
$$\frac{dS_t}{S_t} = r\,dt + \sqrt{v_t}\,dW_t$$
 
$$dv_t = k(\theta - v_t)\,dt + \eta\sqrt{v_t}\,dB_t$$
 
$$\langle dW_t, dB_t \rangle = \rho\,dt$$
 
where the two correlated Brownian motions are constructed from two independent ones $B_t$ and $Z_t$:
 
$$W_t = \rho B_t + \sqrt{1-\rho^2} Z_t$$
 
| Parameter | Meaning |
|-----------|---------|
| $r$ | Risk-free interest rate |
| $k$ | Speed of mean-reversion of variance |
| $\theta$ | Long-term variance (asymptotic) |
| $\eta$ | Volatility of volatility |
| $\rho$ | Correlation between asset and variance Brownians |
| $v_0$ | Initial variance |
 
The long-term mean of the variance satisfies $\lim_{t \to \infty} \mathbb{E}[v_t] = \theta$.
 
---
 
## Milstein Discretization
 
The asset price and variance are discretized using the **Milstein scheme**:
 
$$S_{i+1} = S_i \exp\left[\left(r - \frac{v_i}{2}\right)\Delta t + \sqrt{v_i}\left(\rho\sqrt{\Delta t}N_1^{(i)} + \sqrt{1-\rho^2}\sqrt{\Delta t}N_2^{(i)}\right)\right]$$
 
$$v_{i+1} = v_i + k(\theta - v_i)\Delta t + \eta\sqrt{v_i}\sqrt{\Delta t}N_1^{(i)} + \frac{\eta^2}{4}\Delta t\left((N_1^{(i)})^2 - 1\right)$$
 
where $N_1, N_2 \stackrel{i.i.d.}{\sim} \mathcal{N}(0,1)$.
 
---
 
# Implemented Functions
 
### Simulation
- `simulate_heston(S0, v0, r, k, theta, eta, rho, T, N, Nmc)` — Simulates $N_{mc}$ trajectories of the asset price and variance using the Milstein scheme. Returns the full paths $\{S_i\}$ and $\{v_i\}$.
---
 
### Option Pricing
- `heston_call_MC(S0, K, v0, r, k, theta, eta, rho, T, N, Nmc)` — Estimates the European Call price using **Estimator 1**:
$$V^{heston}_{est1} = \frac{1}{N_{mc}} \sum_{n=1}^{N_{mc}} \max(S_T^{(n)} - K, 0)$$
- `heston_call_MC_antithetic(S0, K, v0, r, k, theta, eta, rho, T, N, Nmc)` — Estimates the Call price using **Estimator 2** (antithetic variance reduction):
$$V^{heston}_{est2} = \frac{1}{2N_{mc}} \sum_{n=1}^{N_{mc}} \left[\max(S_T^{(n)} - K, 0) + \max(S_T^{sym,(n)} - K, 0)\right]$$
where symmetric trajectories use $(-N_1, -N_2)$ instead of $(N_1, N_2)$.
 
---
 
### Greeks
- `greek_theta_heston(S0, K, v0, r, k, theta, eta, rho, T, N, Nmc, h1)` — Computes the Greek $\partial V / \partial \theta$ by centered finite difference using the **same random numbers**:
$$\frac{\partial V^{heston}}{\partial \theta} = \frac{V^{heston}_{est2}(\theta + h_1) - V^{heston}_{est2}(\theta - h_1)}{2h_1}$$
- `greek_eta_heston(S0, K, v0, r, k, theta, eta, rho, T, N, Nmc, h2)` — Computes the Greek $\partial V / \partial \eta$ by centered finite difference using the **same random numbers**.
---
 
### Black-Scholes (for Newton algorithm)
- `BS_call(S0, K, r, T, sigma)` — Computes the Black-Scholes Call price.
- `vega_BS(S0, K, r, T, sigma)` — Computes the Black-Scholes Vega:
$$\frac{\partial V^{BS}}{\partial \sigma} = S_0 \sqrt{T} \frac{1}{\sqrt{2\pi}} e^{-d_1^2/2}$$
---
 
### Implied Volatility
- `newton_implied_vol(V_heston, S0, K, r, T)` — Finds implied volatility $\sigma^{impl}$ such that $V^{BS}(\sigma^{impl}) = V^{heston}$ using the Newton algorithm:
$$\sigma_{n+1} = \sigma_n - \frac{V^{BS}(\sigma_n) - V^{heston}}{\partial V^{BS}/\partial \sigma(\sigma_n)}$$
with starting point $\sigma_0 = \sqrt{2\left|\frac{\ln(S_0/K) + rT}{T}\right|}$.
---
 
### Calibration
- `levenberg_marquardt_heston(V_market, K_market, T, S0, r, k, v0, rho)` — Calibrates $(\theta, \eta)$ by minimizing:
$$\Phi(\theta, \eta) = \sum_{p=1}^{P} \left|V^{marche}(K_p, T_p) - V^{heston}(K_p, T_p, \theta, \eta)\right|^2$$
---
 
# Part I — Simulation and Log-Return Analysis
 
## Simulation Parameters
 
| Parameter | Value |
|-----------|-------|
| $N_{mc}$ | 10 000 |
| $S_0$ | 1 |
| $K$ | 1 |
| $T$ | 0.5 |
| $r$ | 0.01 |
| $k$ | 2 |
| $\theta$ | 0.04 |
| $\eta$ | 0.3 |
| $v_0$ | 0.04 |
| $\Delta t$ | $T/100$ |
 
## Tasks
 
- Plot 4–5 trajectories of the variance $v_t$
- Plot 4–5 trajectories of the asset price $S_t$
- Simulate the log-return $R = \ln(S_T / S_0)$ for $\rho \in \{0, 0.9, -0.9\}$ and plot its density
**Interpretation of $\rho$:**
- $\rho = 0$ → symmetric distribution, close to normal
- $\rho > 0$ → right heavy tail (volatility increases with price)
- $\rho < 0$ → left heavy tail (volatility increases when price falls) — **typical for equity markets**
---
 
# Part II — Monte-Carlo Call Pricing and Variance Reduction
 
## Estimator 1 vs Estimator 2
 
Both estimators are compared on the same plot for $N_{mc} = 100$ to illustrate the **variance reduction** achieved by the antithetic method.
 
The Call price curve $S_0 \mapsto V(0, S_0)$ is plotted for $S_0 \in [0, 20]$.
 
---
 
# Part III — Greeks by Monte-Carlo
 
## Parameters
 
| Parameter | Value |
|-----------|-------|
| $N_{mc}$ | 1 000 – 10 000 |
| $h_1 = h_2$ | 0.1 |
| $S_0$ | 10 |
| $K$ | $[0, 20]$ |
| $T$ | 0.5 |
| $r$ | 0.1 |
| $k$ | 3 |
| $\rho$ | 0.5 |
| $\theta$ | 0.2 |
| $\eta$ | 0.5 |
| $v_0$ | 0.04 |
 
The curves $K \mapsto \partial V / \partial \theta$ and $K \mapsto \partial V / \partial \eta$ are plotted for $K \in [0, 20]$.
 
---
 
# Part IV — Model Calibration
 
## Fixed Parameters
 
$$k = 3, \quad v_0 = 0.04, \quad \rho = 0.5, \quad S_0 = 10, \quad T = 0.5, \quad r = 0.01$$
 
## Calibrated Parameters
 
$$\beta_1 = \theta, \quad \beta_2 = \eta$$
 
## Market Data
 
| $K_p$ | 8 | 8.4 | 8.8 | 9.2 | 9.6 | 10 | 10.4 | 10.8 | 11.2 | 11.6 |
|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| $V^{marche}_p$ | 2.0944 | 1.7488 | 1.4266 | 1.1456 | 0.8919 | 0.7068 | 0.5461 | 0.4187 | 0.3166 | 0.2425 |
 
| $K_p$ | 12 | 12.4 | 12.8 | 13.2 | 13.6 | 14 | 14.4 | 14.8 | 15.2 | 15.6 | 16 |
|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| $V^{marche}_p$ | 0.1860 | 0.1370 | 0.0967 | 0.0715 | 0.0547 | 0.0381 | 0.0306 | 0.0239 | 0.0163 | 0.0139 | 0.086 |
 
## Algorithm: Levenberg-Marquardt
 
**Initialization:** $\theta^{(0)} = 0.2$, $\eta^{(0)} = 0.5$, $\varepsilon = 10^{-4}$, $\lambda = 0.01$
 
**Update rule:**
 
$$M = J^T J + \lambda I, \qquad d^{(k)} = -M^{-1} J^T r^{(k)}, \qquad \beta^{(k+1)} = \beta^{(k)} + d^{(k)}$$
 
**Stopping criterion:** $\sqrt{\sum_i (d_i^{(k)})^2} < \varepsilon$
 
**Constraints (projection):**
- If $\theta > 1$ or $\theta < 0$ → reset $\theta = \theta^{(0)}$
- If $\eta > 1$ or $\eta < 0$ → reset $\eta = \eta^{(0)}$
---
 
# Part V — Volatility Smile
 
For each strike $K_i \in [5, 5.5, \ldots, 20]$, the implied volatility $\sigma^{impl}_i$ is computed such that:
 
$$V^{BS}(S_0, T, K_i, r, \sigma^{impl}_i) = V^{heston}_i(T, K_i)$$
 
**Arbitrage check** before Newton: the Heston price must satisfy:
 
$$\max(S_0 - Ke^{-rT}, 0) < V^{heston} < S_0$$
 
## Parameters (Standard)
 
| Parameter | Value |
|-----------|-------|
| $S_0$ | 10 |
| $T$ | 0.5 |
| $r$ | 0.1 |
| $k$ | 0.3 |
| $\rho$ | 0.7 |
| $\theta$ | 0.3 |
| $\eta$ | 0.4 |
| $v_0$ | 0.03 |
 
## Parameters (Gatheral)
 
| Parameter | Value |
|-----------|-------|
| $T$ | 0.5 |
| $r$ | 0.1 |
| $k$ | 1.32 |
| $\rho$ | -0.7 |
| $\theta$ | 0.04 |
| $\eta$ | 0.4 |
| $v_0$ | 0.02 |
 
*(Parameters from J. Gatheral, "The Volatility Surface")*
 
---
 
# Libraries Used
 
```
numpy
matplotlib
scipy
```
 
Install dependencies with:
 
```bash
pip install numpy matplotlib scipy
```
 
---
 
# Usage
 
Run the main script:
 
```bash
python TP5.py
```
 
The script will:
1. Simulate Heston trajectories and log-return densities for different $\rho$ (Part I)
2. Price European Calls with Estimator 1 and 2, compare variance reduction (Part II)
3. Compute and plot Greeks $\partial V/\partial \theta$ and $\partial V/\partial \eta$ (Part III)
4. Calibrate $(\theta, \eta)$ via Levenberg-Marquardt on market data (Part IV)
5. Reconstruct the implied volatility smile via Newton algorithm (Part V)
---
 
# Key Concepts Covered
 
- Heston stochastic volatility model
- Milstein discretization scheme
- Monte-Carlo simulation of correlated Brownian motions
- Antithetic variance reduction
- Monte-Carlo Greeks by finite difference (pathwise)
- Levenberg-Marquardt nonlinear least squares calibration
- Implied volatility inversion via Newton algorithm
- Volatility smile and skew interpretation
---
 
# Disclaimer
 
This repository contains the **code and methodology** used in an academic laboratory session at CY Tech (Université Paris-Cergy). The original written report has been removed due to **academic and data usage restrictions**.
 
