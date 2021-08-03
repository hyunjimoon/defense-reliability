import numpy
import matplotlib
import matplotlib.pyplot as plot

def plot_gp_realizations(data, samples):
    """Plot posterior Gaussian process realizations"""

    mid_teal="#487575"

    x = data['x']
    y = data['y']

    N_predict = data['N_predict']
    x_predict = data['x_predict']
    y_predict = data['y_predict']

    cmap = matplotlib.cm.get_cmap('Reds')
    I = len(samples['f_predict'])
    buff = 0.25 * I

    for i in range(I):
        f_predict = samples['f_predict'][i]
        plot.plot(x_predict, f_predict, color = cmap((i + buff) / (I + 2 * buff)), alpha=0.1)

    plot.scatter(x_predict, y_predict, color="white", s=8, zorder=3)
    plot.scatter(x_predict, y_predict, color=mid_teal, s=5, zorder=4)
    plot.scatter(x, y, color="white", s=35, zorder=3)
    plot.scatter(x, y, color="black", s=25, zorder=4)

    plot.gca().set_title("Posterior Realizations")
    plot.gca().set_xlim([-10, 10])
    plot.gca().set_xlabel("x")
    plot.gca().set_ylim([-5, 10])
    plot.gca().set_ylabel("y")

    plot.show()

def plot_gp_low_sigma_gp_realizations(data, samples, sigma_low=0.5):
    """Plot posterior Gaussian process realizations"""

    mid_teal="#487575"

    x = data['x']
    y = data['y']

    N_predict = data['N_predict']
    x_predict = data['x_predict']
    y_predict = data['y_predict']

    cmap = matplotlib.cm.get_cmap('Reds')
    I = len(samples['f_predict'])
    buff = 0.25 * I

    for i in range(I):
        f_predict = samples['f_predict'][i]
        plot.plot(x_predict, f_predict, color="#CCCCCC")

    for i in range(I):
        f_predict = samples['f_predict'][i]
        if samples['sigma'][i] < sigma_low:
            plot.plot(x_predict, f_predict, color = cmap((i + buff) / (I + 2 * buff)), alpha=0.1)

    plot.scatter(x_predict, y_predict, color="white", s=8, zorder=3)
    plot.scatter(x_predict, y_predict, color=mid_teal, s=5, zorder=4)
    plot.scatter(x, y, color="white", s=35, zorder=3)
    plot.scatter(x, y, color="black", s=25, zorder=4)

    plot.gca().set_title("Posterior Realizations (sigma < " + str(sigma_low) + ")")
    plot.gca().set_xlim([-10, 10])
    plot.gca().set_xlabel("x")
    plot.gca().set_ylim([-5, 10])
    plot.gca().set_ylabel("y")

    plot.show()

def plot_gp_post_quantiles(data, samples):
    """Plot posterior Gaussian process quantiles"""
    import matplotlib
    import matplotlib.pyplot as plot

    light="#DCBCBC"
    light_highlight="#C79999"
    mid="#B97C7C"
    mid_highlight="#A25050"
    dark="#8F2727"
    dark_highlight="#7C0000"
    mid_teal="#487575"

    x = data['x']
    y = data['y']

    N_predict = data['N_predict']
    x_predict = data['x_predict']
    y_predict = data['y_predict']

    P = [10, 20, 30, 40, 50, 60, 70, 80, 90]    
    creds = [numpy.percentile(samples['f_predict'][:,n], P) for n in range(N_predict)]
                
    plot.fill_between(x_predict, [c[0] for c in creds], [c[8] for c in creds],
                      facecolor=light, color=light)
    plot.fill_between(x_predict, [c[1] for c in creds], [c[7] for c in creds],
                      facecolor=light_highlight, color=light_highlight)
    plot.fill_between(x_predict, [c[2] for c in creds], [c[6] for c in creds],
                      facecolor=mid, color=mid)
    plot.fill_between(x_predict, [c[3] for c in creds], [c[5] for c in creds],
                     facecolor=mid_highlight, color=mid_highlight)
    plot.plot(x_predict, [c[4] for c in creds], color=dark)

    plot.scatter(x_predict, y_predict, color="white", s=8, zorder=3)
    plot.scatter(x_predict, y_predict, color=mid_teal, s=5, zorder=4)
    plot.scatter(x, y, color="white", s=35, zorder=3)
    plot.scatter(x, y, color="black", s=25, zorder=4)
 
    plot.gca().set_title("Posterior Quantiles")
    plot.gca().set_xlim([-10, 10])
    plot.gca().set_xlabel("x")
    plot.gca().set_ylim([-5, 10])
    plot.gca().set_ylabel("y")

    plot.show()

def plot_gp_retro_pred_quantiles(data, samples):
    """Plot posterior retrodictive Gaussian process quantiles"""
    import matplotlib
    import matplotlib.pyplot as plot

    light="#DCBCBC"
    light_highlight="#C79999"
    mid="#B97C7C"
    mid_highlight="#A25050"
    dark="#8F2727"
    dark_highlight="#7C0000"
    mid_teal="#487575"

    x = data['x']
    y = data['y']

    N_predict = data['N_predict']
    x_predict = data['x_predict']
    y_predict = data['y_predict']

    P = [10, 20, 30, 40, 50, 60, 70, 80, 90]    
    creds = [numpy.percentile(samples['y_predict'][:,n], P) for n in range(N_predict)]
                
    plot.fill_between(x_predict, [c[0] for c in creds], [c[8] for c in creds],
                      facecolor=light, color=light)
    plot.fill_between(x_predict, [c[1] for c in creds], [c[7] for c in creds],
                      facecolor=light_highlight, color=light_highlight)
    plot.fill_between(x_predict, [c[2] for c in creds], [c[6] for c in creds],
                      facecolor=mid, color=mid)
    plot.fill_between(x_predict, [c[3] for c in creds], [c[5] for c in creds],
                     facecolor=mid_highlight, color=mid_highlight)
    plot.plot(x_predict, [c[4] for c in creds], color=dark)

    plot.scatter(x_predict, y_predict, color="white", s=8, zorder=3)
    plot.scatter(x_predict, y_predict, color=mid_teal, s=5, zorder=4)
    plot.scatter(x, y, color="white", s=35, zorder=3)
    plot.scatter(x, y, color="black", s=25, zorder=4)

    plot.gca().set_title("Posterior Retrodictive Quantiles")
    plot.gca().set_xlim([-10, 10])
    plot.gca().set_xlabel("x")
    plot.gca().set_ylim([-5, 10])
    plot.gca().set_ylabel("y")

    plot.show()
