# Reliabile Defense System
[NextOpt team](https://www.hyunjimoon.com/blog/vision-of-my-startup-nextopt/) build a reliable defense system with failure prediction, system risk measure and optimization, mechanism design techniques.

## Target 
### 1. Posterior of Quantity of interest 
Quantity of interest (phi) is similar to hyperparameter as it affects data only by means of model parameter (theta) and also they both await modeler's decision. Retroactive exploration of data generation process phi -> theta -> (x,y) process (e.g. decision -> parameter -> data) provides tools for (X, Y) -> Theta -> Phi and therefore, posterior of Phi. Captitalized (X, Y, Theta) represent population random variable (with identified distribution) while (x, y, theta) denote sample values. Two approaches exist for phi posterior:
- Deductive `Phi|X,Y` density via analytic calculation i.e. get `Phi|Theta` and `Theta|X,Y` then marginalizing out `Theta`.  
- Inductive `phi|x,y` sample via computation. Satisfying targeted simulation-based calibration is the necessary condition which determines sample credibility by inspecting the consistency of sampling mechanism. Sampling mechanism consists of three simulators: prior, data, and posterior.

![image](https://user-images.githubusercontent.com/30194633/128450710-4b41ff94-0026-4ff0-b037-db6196d05a7b.png)

### 2. System risk and interaction 
- Systme risk management on top of its **vertical and horizontal interaction** is our interest and we use **hierarchical and mixture model** as our main frame.
- System dynamics simulation approach is being updated in Moon S. [vensim](https://m.blog.naver.com/mseongam/222059202785) blog.

### 3. Mechanism Design
- Contract is being updated in Moon S. [supply contract](https://m.blog.naver.com/PostList.naver?blogId=mseongam&categoryNo=6&listStyle=style1) blog.

--- 
## Tool 
### 1. Y|y, X|x, Theta|XY, Phi|Theta
In `R` folder, topics with the following keyword are sorted. Notice the change of background space: from **data = Y** to **parameter = Theta** than to **quantity of interest := f(Theta) | X**. X is a predictor with (n x p) matrix and y is a data with length n vector. Decision is the best example for quantity of interest (QI)  whose posterior distribution is attained with the help of parameter.

1. `Y|y` in data space. 
In the presence of limited data, impute raw data with certified assumptions and construct the generative process from Theta to E[Y].
 - transform scale and distribution 
 - explore data via conditioning predictor: E[theta|X=a] vs E[theta|X=b]
 - set resolution. Bin the time axis (e.g. age of the product) by considering the amount of data for individual period interval and the granuality of required forecast
 - identify outlier. Drop noise that murk the main relationship between Ey and theta while preserve the main component of thier relationship Impute two types of data. Data exploration and signal to noise ratio prior knowledge are helpful. Quantile-based drop (or replacement) is the most common.
 - incorporate expert knowledge on distribution and range of data. 
 - construct test set for each scenario or layer. For K-layer hierarchical model (HM), K possible cases exist that needs separate testing; e.g two-layer HM with 1-5-99 engine_archetype(phi)-engine(theta[1..5])-ship(mu[1..99]). Two testsets, first with known engine and unkown ship, and the second, both unkown engine and ship, need construction. 
 
2. `X|x` in parameter space.
 - generate scaled timeseries features: trend, seasonality, event, self-lag etc
 - generate hierarchical feature i.e. group index 
 - select feature. e.g. blackbox forward and backward selection algorithm while more adaptive spike-and-slad or more transparent causal effect based selection are possible                                                  

3. `Theta|XY` is from data to parameter space. 
 - infer parameter values given data for each predictor
 - design pooling structure between different predictors with the assumption Theta|X = a is similar to Theta|X = b
 - `extreme` increase estimation/simulation efficiency using splitting, exploiting regeneration structure**, verification techniques to model extreme event where `Theta|XY` is highly inefficient 

4. `Phi|Theta` is from parameter to quantity of interest space. 
 - marginalize out nuisance parameter to calculate posterior of continuous QI (e.g. scaled failure counts) or marginal likelihood of discrete QI (e.g. preventive maintenance period)s: schedule two types of maintenance: preventive triggered by inspection & corrective 
 - remove the middle target (prediction) to directly address the decision problem ([Bayesian optimization](https://ieeexplore.ieee.org/document/7352306), [Smart predict then optimize](https://www.ima.umn.edu/materials/2018-2019.1/W10.3-5.18/27490/SPO_121317.pdf))  
 - design pooling strucutre in data space: determine the model weight for Bayesian model averaging and stacking especially in HM
 - design pooling strucutre in parameter space: aggregate parameter distribution from different models on joint parameter space
 - identify QI options based system requirement e.g. QI = preventive maintenance triggered by inspection & corrective 

### 2. Verification and Validation
- To verify `theta|x,y` or `qi|x,y`, [Simulation-based calibration](https://mc-stan.org/docs/2_27/stan-users-guide/simulation-based-calibration.html) which is maintained in another [repo](https://github.com/hyunjimoon/SBC/tree/api-variant) can be applied.

### References:
#### Modeling
- *[Mixed pooling of seasonality for time series forecasting: An application to pallet transport data](https://www.researchgate.net/publication/346259196_Mixed_pooling_of_seasonality_for_time_series_forecasting_An_application_to_pallet_transport_data) Moon, H., Song, B., & Lee H. (2020), under revision.

- **[Exploiting regenerative structure to estimate finite time averages via simulation](http://www.columbia.edu/~ww2040/WanmoRevised.pdf) Kang, W., Shahabuddin, P., & Whitt, W. (2006)

#### Computation 
Rare event estimation and simulation techniques

- [Modelling Extremal Events estimation](https://www.springer.com/gp/book/9783540609315) explains how to estimate the tails of distributions and Ch.6 is the most illustrative.

- [Introduction to Rare Event Simulation](https://www.springer.com/gp/book/9780387200781) shows efficient Monte Carlo computation to estimate occurrence proabaility of rare events. [Rare-Event Simulation Techniques: An Introduction and Recent Advances (survey)](https://www.sciencedirect.com/science/article/pii/S092705070613011X) introduce a shorter overview.

#### System risk
[An Axiomatic Approach to Systemic Risk](https://www.jstor.org/stable/23443854?seq=1#metadata_info_tab_contents)

#### Domain knowledge of (Korean) military 
- [Naval Vessel Spare Parts Demand Forecasting Using Data Mining, Yoon17](http://www.ksie.ne.kr/journal/article.php?code=58051)
- [Forecasting Spare Parts Demand of Military Aircraft: Comparisons of Data Mining Techniques and Managerial Features from the Case of South Korea](https://www.mdpi.com/2071-1050/12/15/6045)
- [모듈형 엔진의 수명관리를 고려한 항공기-임무 할당 모형](https://nextoptext.slack.com/archives/C013D35MN9J/p1610967442002200)
![image](https://user-images.githubusercontent.com/30194633/112784366-79e7fd00-908c-11eb-9aab-3728108675c7.png)

#### Contract and behavior
tbc

##### Note
- For the techniques we (NextOpt team) developed, refer to the writeup folder.
- Inventory management is within our research scope but, not in the near future.
