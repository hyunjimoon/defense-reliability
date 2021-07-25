# Reliabile Defense System
This repository includes modeling approaches and references on building a reliable defense system. Not only the Bayesian workflow methodologies and its modification to system structure but also the application to real-world problems such as part failure prediction, maintenance policy optimization are our contribution. Topics are categorized with `keyword_`. Notice the change of background space: from **data = Y** to **parameter = Theta** than to **QI := f(Theta) | X**. X is a predictor with (n x p) matrix and y is a data with length n vector. QI is a quantity of interest which is the function of inferred parameters; the absense of its closed form can be offset by its evaluation ability. evaluated value of the function at every xi. QI can be offset by its ability to evaluate the value for every point. Captitalized (Y, Theta) represent population random variable (with identified distribution) while the other (y, theta) denote sample values (given data).

1. `Y_bar_y` in data space. 
In the presence of limited data, impute raw data with certified assumptions and construct the generative process from Theta to E[Y].
 - scaling (todo)
 - explore data via conditioning predictor: E[theta|X=a] vs E[theta|X=b]
 - set resolution. Bin the time axis (e.g. age of the product) by considering the amount of data for individual period interval and the granuality of required forecast
 - identify outlier. Drop noise that murk the main relationship between Ey and theta while preserve the main component of thier relationship Impute two types of data. Data exploration and signal to noise ratio prior knowledge are helpful. Quantile-based drop (or replacement) is the most common.
 - incorporate expert knowledge on distribution and range of data 
 - construct test set for each scenario. For hierarchical model (HM), each layer needs separate examination; for HM where information is pooled among the engines that share and among the ships that share the engine among the ships  shared engine layer, two testship with out of sample engine   insample engine but  engine with be examined. out of sample type or known exists but unobserved and doesn't exist (e.g. failure of ). 

2. `X_bar_x` in parameter space. (todo refactoring)
 - generate scaled timeseries features: trend, seasonality, event, self-lag etc.
 - generate hierarchical feature i.e. group index 
 - select feature. e.g. blackbox forward and backward selection algorithm while more adaptive spike-and-slad (todo) or more transparent causal effect based selection (todo) are possible                                                  

3. `Theta_bar_Y_bar_X` is from data to parameter space. 
 - infer parameter values given data for each predictor
 - design pooling structure between different predictors with the assumption: Theta_bar_Y|X = a is similar to Theta_bar_Y|X = b
 - `extreme` increase estimation/simulation efficiency using splitting, exploiting regeneration structure**, verification techniques to model extreme event where `Theta_bar_Y_bar_X` is highly inefficient 

4. `QI_bar_Theta` is from parameter to QI space.
 - marginalize out nuisance parameter to calculate posterior of continuous QI (e.g. scaled failure counts) or marginal likelihood of discrete QI (e.g. preventive maintenance period)s: schedule two types of maintenance: preventive triggered by inspection & corrective 
 - remove the middle target (prediction) to directly address the decision problem ([Bayesian optimization](https://ieeexplore.ieee.org/document/7352306), [Smart predict then optimize](https://www.ima.umn.edu/materials/2018-2019.1/W10.3-5.18/27490/SPO_121317.pdf))  
 - design pooling strucutre in data space: determine the model weight for Bayesian model averaging and stacking especially in HM (todo)
 - design pooling strucutre in parameter space: aggregate parameter distribution from different models on joint parameter space (todo)
 - identify QI options based system requirement e.g. QI = preventive maintenance triggered by inspection & corrective 

### References:
#### Modeling
- *[Mixed pooling of seasonality for time series forecasting: An application to pallet transport data](https://www.researchgate.net/publication/346259196_Mixed_pooling_of_seasonality_for_time_series_forecasting_An_application_to_pallet_transport_data) Moon, H., Song, B., & Lee H. (2020), under revision.

- **[Exploiting regenerative structure to estimate finite time averages via simulation](http://www.columbia.edu/~ww2040/WanmoRevised.pdf) Kang, W., Shahabuddin, P., & Whitt, W. (2006)

#### Computation 
Rare event estimation and simulation techniques

- [Modelling Extremal Events estimation](https://www.springer.com/gp/book/9783540609315) explains how to estimate the tails of distributions and Ch.6 is the most illustrative.

- [Introduction to Rare Event Simulation](https://www.springer.com/gp/book/9780387200781) shows efficient Monte Carlo computation to estimate occurrence proabaility of rare events. [Rare-Event Simulation Techniques: An Introduction and Recent Advances (survey)](https://www.sciencedirect.com/science/article/pii/S092705070613011X) introduce a shorter overview.

#### Domain knowledge of (Korean) military 
- [Naval Vessel Spare Parts Demand Forecasting Using Data Mining, Yoon17](http://www.ksie.ne.kr/journal/article.php?code=58051)
- [Forecasting Spare Parts Demand of Military Aircraft: Comparisons of Data Mining Techniques and Managerial Features from the Case of South Korea](https://www.mdpi.com/2071-1050/12/15/6045)
- [모듈형 엔진의 수명관리를 고려한 항공기-임무 할당 모형](https://nextoptext.slack.com/archives/C013D35MN9J/p1610967442002200)
![image](https://user-images.githubusercontent.com/30194633/112784366-79e7fd00-908c-11eb-9aab-3728108675c7.png)

##### Note
- For the techniques we (NextOpt team) developed, refer to the writeup folder.
- Inventory management is within our research scope but, not in the near future.
