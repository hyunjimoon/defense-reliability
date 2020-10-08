import numpy as np
def inverse_transform(x, lmbda, mean, std, resolution=0):
    """
    sklearn powertransform(yeo-johnson) inverse transformer
    :param lmbda: lambda value of Yeo-Johnson Transform (PowerTransformer.lambdas_)
    :param mean: mean value of StandardScaler (PowerTransformer._scaler.mean_)
    :param std: StandaScaler standard deviation (PowerTransformer._scaler.scale_)
    :param resolution: np.round resolution. Keep at 0 unless necessary.
    ex: 
    lambda = 0.21649144
    mean = 4.96401329
    std = 2.3079645
    """
    
    x = x * std + mean
    
    x_inv = np.zeros_like(x)
    pos = x >= 0

    # when x >= 0
    if abs(lmbda) < np.spacing(1.):
        x_inv[pos] = np.exp(x[pos]) - 1
    else:  # lmbda != 0
        x_inv[pos] = np.power(x[pos] * lmbda + 1, 1 / lmbda) - 1

    # when x < 0
    if abs(lmbda - 2) > np.spacing(1.):
        x_inv[~pos] = 1 - np.power(-(2 - lmbda) * x[~pos] + 1,
                                   1 / (2 - lmbda))
    else:  # lmbda == 2
        x_inv[~pos] = 1 - np.exp(-x[~pos])

    return np.round(x_inv, resolution)