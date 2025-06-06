import yfinance as yf
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import math

start_date = datetime(2023, 6, 8)
end_date = start_date + timedelta(days=366)

##### Download data ######
# First 4 is the benchmark data
tickers =["CBA.AX", "NAB.AX", "WBC.AX", "ANZ.AX","^AXJO"]
benchmarkTicker = "^AXJO"
securities = yf.download(tickers=tickers, start=start_date, end=end_date)

# Calculate returns per ticker
for ticker in tickers:
    securities[("Returns",ticker)] = securities[("Close",ticker)].pct_change()

initial_value = 2500

# Caculate the accumulative value of the individual holdings
for ticker in tickers:
    # If you want to include the starting value (i.e., first row as initial value),
    # then fill NaN return with 0 for first row to start at initial_value
    returns = securities[("Returns", ticker)].fillna(0)

    # Compute value series
    value_series = initial_value * (1 + returns).cumprod()

    # Store it in the same MultiIndex DataFrame
    securities[("Value", ticker)] = value_series


##### Calculate portfolio returns ######

# Extract the portfolio tickers
columns = []
for i in range(len(tickers)-1):
    ticker = tickers[i]
    columns.append(('Value', ticker))

# Create new data frame containing just the return columns
portfolioIndividualValues = securities[columns]
portfolioIndividualValues["Total Value"] = securities[columns].sum(axis=1)
portfolioIndividualValues["Returns"] = portfolioIndividualValues["Total Value"].pct_change()


returns_df = pd.DataFrame({
    "Portfolio Returns": portfolioIndividualValues["Returns"],
    "Benchmark Returns": securities[("Returns",tickers[-1])]
}, index=securities.index)

# Create a differenced series
returns_df["Difference"] = returns_df["Portfolio Returns"] - returns_df["Benchmark Returns"]

# Plot the portfolio returns
returns_df["Portfolio Returns"].plot(xlabel="Date", 
                                     ylabel="Net Daily Return", 
                                     title="Daily Returns of Big 4 Banks Portfolio 1 year from 08-06-2023",
                                     figsize=(10,6),
                                     color="blue")
plt.savefig("Portfolio_returns.png")
#plt.show()

# Plot the benchmark returns
returns_df["Benchmark Returns"].plot(xlabel="Date", 
                                     ylabel="Net Daily Return",
                                     title="Daily Returns of ASX index 1 year from 08-06-2023",
                                     figsize=(10,6),
                                     color="orange")
plt.savefig("Benchmark_returns.png")
#plt.show()

returns_df[["Portfolio Returns","Benchmark Returns"]].plot(xlabel="Date", 
                                     ylabel="Net Daily Return",
                                     title="Daily Returns of ASX index and Portfolio 1 year from 08-06-2023",
                                     figsize=(10,6),
                                     color=["blue","orange"])
plt.savefig("returns.png")
#plt.show()

returns_df["Difference"].plot(xlabel="Date", 
                              ylabel="Net Daily Return Difference",
                              title="Portfolio returns less Benchmark returns from 08-06-2023",
                              figsize=(10,6),
                              color="green")
plt.savefig("Returns_comparison.png")
#plt.show()

##### Calculate various statistics ######

# Porfolio Statistics
meanPortfolio = returns_df["Portfolio Returns"].mean()
variancePortfolio = returns_df["Portfolio Returns"].var()
standardDeviationPortfolio = returns_df["Portfolio Returns"].std()
# variancePortfolio = returns_df.apply(lambda row: (row["Portfolio Returns"]-meanPortfolio)**2,axis=1).mean()
# standardDeviationPortfolio = math.sqrt(variancePortfolio)

annualMeanPortfolio = meanPortfolio*252
annualVariancePortfolio = variancePortfolio * 252
standardDeviationAnnualPortfolio = standardDeviationPortfolio*math.sqrt(252)

# Benchmark Statistics
meanBenchmark = returns_df["Benchmark Returns"].mean()
varianceBenchmark = returns_df["Benchmark Returns"].var()
standardDeviationBenchmark = returns_df["Benchmark Returns"].std()
# varianceBenchmark = returns_df.apply(lambda row: (row["Benchmark Returns"]-meanBenchmark)**2,axis=1).mean()
# standardDeviationBenchmark = math.sqrt(variancePortfolio)

annualMeanBenchmark = meanBenchmark*252
annualVarianceBenchmark = varianceBenchmark * 252
standardDeviationAnnualBenchmark = standardDeviationBenchmark*math.sqrt(252)

# Other statistics
standardDeviationDifferenceAnnual = returns_df["Difference"].std()*math.sqrt(252)

covariance = returns_df.cov().loc["Portfolio Returns", "Benchmark Returns"]
#covariance = returns_df.apply(lambda row: (row["Portfolio Returns"]-meanPortfolio)*
#                              (row["Benchmark Returns"]-meanBenchmark),axis=1).mean()
annualCovariance = covariance*252


# Portfolio performance metrics

beta = annualCovariance / (annualVarianceBenchmark)

riskFreeRate = 0.03

jenAlpha = annualMeanPortfolio - (riskFreeRate+beta*(annualMeanBenchmark-riskFreeRate)) 

sharpeRatio = (annualMeanPortfolio-riskFreeRate) / standardDeviationAnnualPortfolio

informationRatio = (annualMeanPortfolio-annualMeanBenchmark) / (standardDeviationDifferenceAnnual)

print(f"Portfolio:\n expected annualised return: {annualMeanPortfolio}\n annualised standard deviation: {standardDeviationAnnualPortfolio}")
print(f"Benchmark:\n expected annualised return: {annualMeanBenchmark}\n annualised standard deviation: {standardDeviationAnnualBenchmark}")
print(f"Performance metrics\n Beta: {beta}\n Jenson's alpha: {jenAlpha}\n Sharpe Ration: {sharpeRatio}\n Information Ratio: {informationRatio}")
