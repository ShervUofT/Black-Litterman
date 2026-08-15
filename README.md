
# Black-Litterman Portfolio Optimization

## Overview

This project explores how the **Black-Litterman portfolio model** can be used to build and evaluate investment portfolios using industry-level market data.

The goal was to compare Black-Litterman portfolios with more traditional portfolio construction methods and to test how different types of investment views affect portfolio performance.

## What the Project Does

The project uses historical Fama-French industry data to:

* Build Black-Litterman portfolios
* Compare them with traditional mean-variance portfolios
* Test different methods for forming investment views
* Backtest the strategies over time
* Evaluate risk, return, diversification, and portfolio stability
* Account for portfolio turnover and transaction costs

## View Strategies

Three different approaches are used to create investment views:

### Trailing Returns

Uses recent relative performance between selected industries to form portfolio views.

### Momentum

Assumes industries that have recently outperformed similar industries may continue to outperform.

### Mean Reversion

Assumes unusually large differences in performance between industries may eventually reverse.

## Portfolio Comparison

The analysis compares several approaches:

* Black-Litterman portfolios
* Traditional mean-variance portfolios
* Portfolios based primarily on investment views
* A market-size weighted benchmark

Performance is evaluated using measures such as return, volatility, Sharpe ratio, portfolio concentration, turnover, and cumulative performance.

## Files

* `CORE_BL_TRAILING.R` — main Black-Litterman model using trailing-return views and the full portfolio comparison
* `BL_MOMENTUM.R` — Black-Litterman model using momentum-based views
* `BL_MEANREVERSION.R` — Black-Litterman model using mean-reversion views

## Data

The project uses historical industry portfolio and market data from the **Kenneth French Data Library**.

## Tools

R, portfolio optimization, Black-Litterman modeling, backtesting, quantitative finance
