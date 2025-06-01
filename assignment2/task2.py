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

# Create lagged close and calculate returns
for ticker in tickers:
    securities[("Close_lag1",ticker)] = securities[("Close",ticker)].shift(1)
    securities[("Returns",ticker)] = (securities[("Close",ticker)] - securities[("Close_lag1",ticker)]) / securities[("Close_lag1",ticker)]

    # Plots    
    # securities[("Returns",ticker)].plot(title=f"{ticker} - 1 Year Returns")
    # plt.show()
    # securities[("Returns",ticker)].plot.hist(bins=100,title=f"{ticker} - 1 Year Returns")
    # plt.show()

# Drop the first day (wont have a return data for the first day)
securities.dropna(inplace=True)

##### Calculate portfolio returns ######

# Extract the portfolio tickers
columns = []
for i in range(len(tickers)-1):
    ticker = tickers[i]
    columns.append(('Returns', ticker))

# Create new data frame containing just the return columns
portfolioIndividualReturns = securities[columns]

# Calculate the portfolio returns
weights = np.array([0.25]*4)
returns_matrix = portfolioIndividualReturns.to_numpy()
returns_matrix.transpose()

# Create a new data frame including the portfolio returns, benchmark returns and difference of the returns
# Alternatively portfolioIndividualReturns.apply(lambda row: row.mean(),axis=1).mean()
portfolioReturns = weights @ returns_matrix.T
returns_df = pd.DataFrame({
    "Portfolio Returns": portfolioReturns,
    "Benchmark Returns": securities[("Returns",tickers[-1])]
}, index=securities.index)

# Create a differenced series
returns_df["Difference"] = returns_df["Portfolio Returns"] - returns_df["Benchmark Returns"]

# Plot the portfolio returns
returns_df["Portfolio Returns"].plot(xlabel="Date", 
                                     ylabel="Net Daily Return", 
                                     title="Daily Returns of Big 4 Banks Portfolio 1 year from 08-06-2023",
                                     figsize=(10,6))

plt.show()

# Plot the benchmark returns
returns_df["Benchmark Returns"].plot(xlabel="Date", 
                                     ylabel="Net Daily Return",
                                     title="Daily Returns of ASX index 1 year from 08-06-2023",
                                     figsize=(10,6))
plt.show()


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
standardDeviationAnnualPortfolio = annualVarianceBenchmark*math.sqrt(252)


# Other statistics
standardDeviationDifferenceAnnual = returns_df["Difference"].std()*math.sqrt(252)

covariance = returns_df.cov().loc["Portfolio Returns", "Benchmark Returns"]
#covariance = returns_df.apply(lambda row: (row["Portfolio Returns"]-meanPortfolio)*
#                              (row["Benchmark Returns"]-meanBenchmark),axis=1).mean()
annualCovariance = covariance*252


# Portfolio performance metrics

beta = annualCovariance / (variancePortfolio*252)

riskFreeRate = 0.03

jenAlpha = annualMeanPortfolio - (riskFreeRate+beta*(annualMeanBenchmark-riskFreeRate)) 

sharpeRatio = (annualMeanPortfolio-riskFreeRate) / standardDeviationAnnualPortfolio

informationRatio = (annualMeanPortfolio-annualMeanBenchmark) / (standardDeviationDifferenceAnnual)

print(f"Beta : {beta}\nJenson's alpah {jenAlpha}\nSharpe Ration {sharpeRatio}\nInformation Ratio {informationRatio}")
