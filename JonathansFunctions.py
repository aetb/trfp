import numpy as np
import pandas as pd


def plt_unix_time_to_CST(ax):
    """I took this function from Alec, it converts the time axis into "real time"
    """
    ax.locator_params(axis='x', nbins=5)
    xticks = ax.get_xticks()
    ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                        for tm in xticks], rotation=0, fontdict={'size':10, 'family':'serif'})

def difference(y, avgt):
    """Bins an array into a specified bin size. ***DIFFERENCE IS A MISLEADING NAME***
    
    This function is not used in the final analysis. It is only used to get an idea of
    what the data looks like. In the final analysis, MovingAvg is used. The last bin 
    if the bin size does not divide the length of the aray, then the last bin (smaller)
    will not be in the binned array (See example)
    
    Args:
        y: the array to bin
        avgt: The bin size
        
    Returns:
        t: an array containing the same number of elements as avgdiff that bookkeeps
           the indices averaged.
        avgdiff: The binnned array.
        example:
        
        difference([0,2,4,6,8], 3)
        outputs
        t = np.array([3])
        avgdiff = np.array([2])
    """
    t = np.arange(int(avgt) , len(y), int(avgt))
    N = len(t)
    Sum = 0
    avgdiff = np.array([])
    for i in range(len(y)):
        Sum += y[i]
        if i % avgt == 0 and i > 0:
            avg = Sum/avgt
            Sum = 0
            avgdiff = np.append(avgdiff, avg)
            
        
    return (t,avgdiff)

def MovingAvg(x, bin_size):
    """Averages data with overlapping bins.
    *** UPDATE: I NOW DO NOT USE THIS FUNCTION AS IT IS SLOW. I NOW USE A PANDAS
    FUNCTION ***
    This function is used in the final data analysis. This function is quite slow.
    Alec has a function that does this much quicker (I believe within Pandas).
    Note that less data is "lost" compared to the difference function.
    
    
    Args:
    x: The array to bin
    bin_size: The overlapping bin size.
    
    Returns:
    avg: the binned array
    example:
    
    MovingAvg([0,1,2,3,4], 2)
    
    returns
    np.array([.5, 1.5, 2.5, 3.5])
    """
    avg = np.array([])
    for i in np.arange(len(x) - bin_size + 1):
    #while (i <= len(x) - bin_size):
        Sum = 0
        j = 0
        for j in np.arange(bin_size):
        #while (j < bin_size):
            Sum += x[i + j]
            #j += 1
        avg = np.append(avg, Sum)
        if(i % 86400 == 0):
            print(str(i) + "iterations done")
        i +=1
    return avg/bin_size



def line(x,m,b): #just a line equation
    return m*x + b

def Multiplot(x_data, y_data, height, width, dimension, sublabels = None, x_axis_labels = None, y_axis_labels = None, figname = None, CST = False, daylines = False, day_offset = 0, logscale = False, day_label = ""):
    """ Makes a plot with given data.
    This function has hard syntax and imputs. Generally not used unless if in a 
    zoom meeting. I used this before I became more comfortable with looping through
    rank n arrays. This is here just because I did not want to erase it. I might delete it 
    later
    """
    #x_data is an array containing arrays with the different data to plot
    #y_data ""
    #width is the width of the whole figure in inches
    #height is the height of the whole figure in inches
    #dimension is an array containing the dimensions of the plot figure np.array([3,2])  would correspond to 6 plots (3 rows, 2 cols)
    #x_axis_labels labels is an array containing the labels on the subplots
    #Similarly for y
    #sublabels are the titles are the subplots (1 x n array of strings), can just be 1 x 1 if all the same
    #Figname is the name of the plot
    #CST will convert the x-axis data into CST time
    #Daylines will plot red lines to show days.
    #day_offset in case we want the 24 hour mark to be somewhere else other than just 24 hours from the start
    for i in range(len(x_data)):
        if len(x_data[i]) != len(y_data[i]):
            print("The dimensions of x[" + str(i) + "] and y[" + str(i) + "] do not match")
    
    equal_x_labels = True
    for k in range(len(x_axis_labels)):
        for l in range(len(x_axis_labels)):
            if k != l:
                if  x_axis_labels[k] != x_axis_labels:
                    equal_x_labels = False
                    
    equal_y_labels = True
    for k in range(len(y_axis_labels)):
        for l in range(len(y_axis_labels)):
            if k != l:
                if  y_axis_labels[k] != y_axis_labels:
                    equal_y_labels = False
    day = np.array([86400,86400])
    fig, axs = plt.subplots(dimension[0], dimension[1], constrained_layout=True)
    fig.set_figheight(height)
    fig.set_figwidth(width)
    for i in range(dimension[0]):
        for j in range(dimension[1]):
            axs[i,j].plot(x_data[dimension[1] * i + j], y_data[dimension[1] * i + j])
            if equal_x_labels:
                axs[i,j].set_xlabel(x_axis_labels[0])
            else:
                axs[i,j].set_xlabel(x_axis_labels[dimension[1] * i + j])
            if equal_y_labels:
                axs[i,j].set_ylabel(y_axis_labels[0])
            else:
                axs[i,j].set_ylabel(y_axis_labels[dimension[1] * i + j])
            axs[i,j].set_title(sublabels[dimension[1] * i + j])
            if daylines:
                days = int((x_data[dimension[1] * i + j].max() - x_data[dimension[1] * i + j].min()) / 86400)
                for n in range(days):
                    axs[i,j].plot(day * (n+1) + day_offset + x_data[dimension[1] * i + j][0], [y_data[dimension[1] * i + j].min(), y_data[dimension[1] * i + j].max()], 'r', label = day_label)
                    axs[i,j].legend()
            if logscale:
                axs[i,j].set_xscale('log')
                axs[i,j].set_yscale('log')
            if CST:
                plt_unix_time_to_CST(axs[i,j])
    fig.suptitle(figname)

def ArrayRound(array, value):
    """ Returns the element closest to the one provided in a non-decreasing array of values.
    
    This function may be changed later to return the index instead of the value.
    NOTE: in case of two values being equally close, the element with lower index
          (i.e. the first one) will be returned.
    
    Args:
        array: array with non-decreasing values (mostly used for arrays to be plotted,
        such as x-axis)
        value: the value that will be searched for.
    Returns:
        An element of the array closest to the value provided.
        example:
        
        ArrayRound([1,2,3,4], 2.5)
        
        returns
        2
    """
    inf = -1
    sup = -1
    if array[0] >= value:
        return array[0]
    else:
        for number in array:
            if number <= value:
                inf = number
        if inf == value:
            return value
        else:
            sup = array[np.where(array == inf)[0] +1]
            if abs(inf - value) < abs(sup - value):
                return inf
            elif abs(inf-value) > abs(sup - value):
                return sup
            else:
                print("The closest the value provided is exactly in between two elements of the array, the lower element has been returned.")
                return inf
        
    
def DayLineFit(time_data, y_data, start):
    """Fits lines through every pair of points that are 24 hours apart beginning at some 
    starting value.
    
    Note that thiss uses pairs of points and not three or more points.
    Also, the start is usually time_data[0]. An ending point can also be added
    in future versions.
    This function is useful to study diurnal cycles/remove drifts.
    
    Args:
        time_data: An array of times in SECONDS. The values of this array must be spaced
                   equally.
        y-data: The data to which lines shall be fitted.
    Returns:
        slopes: the set of all slopes obtained from the fits, with the first element being 
                the fit where the start value is fit
        intercepts: similar to slopes but with intercepts.
        
        Note: these arrays are most likely going to be averaged.
        example: try using this on sin(2*pi*t/86400) + t + 4.
                 Should return an array of slopes and
                intercepts that average to 1 and 4 respectively
                (or faily close to those values)
    """
    step = ArrayRound(time_data, start)
    difference = time_data[1] - time_data[0]
    day = int(round(86400/difference,0))
    slopes = np.array([])
    intercepts = np.array([])
    test = np.array([step, step + (day*difference)])
    while step + 86400 <= time_data[-1]:
        x = np.array([step, step + (day*difference)])
        y = np.array([y_data[np.where(time_data == step)[0]], y_data[np.where(time_data == step + (day*difference))[0]]])
        z = np.polyfit(x,y,1)
        slopes = np.append(slopes, z[0])
        intercepts = np.append(intercepts, z[1])
        step = time_data[np.where(time_data == step)[0][0] +1]
    return slopes, intercepts

def DayParFit(time_data, y_data, start):
    """This function is identical to DayLineFit, except parabolic.
    
    returns: ax^2 + bx + c
        a: second-degree coefficients
        b: first-degree coefficients
        c: y-intercepts
    """
    #Fits parabola on points 24 hours apart. Assuming the time data is in seconds, takes into account any binning done
    #Assumes time_data points are equally spaced; hence the the difference variable
    step = ArrayRound(time_data, start)
    difference = time_data[1] - time_data[0]
    day = int(round(86400/difference,0))
    a = np.array([])
    b = np.array([])
    c = np.array([])
    while step + 2*86400 <= time_data[-1]:
        x = np.array([step, step + (day*difference), step + (2*(day*difference))])
        y = np.array([y_data[np.where(time_data == step)[0]], y_data[np.where(time_data == step + (day*difference))[0]], y_data[np.where(time_data == step + (2*(day*difference)))[0]]])
        z = np.polyfit(x,y,2)
        a = np.append(a, z[0])
        b = np.append(b, z[1])
        c = np.append(c, z[2])
        step = time_data[np.where(time_data == step)[0][0] +1]
    return a, b, c

def LogSlope(x,y,a): #Returns the slope of a graph at x = a, a is the actual number, not log of the point of interest.
    figwidth = 13
    figheight = 5
    a = ArrayRound(x,a)
    t = np.where(x == a)[0][0]
    diff_y32 = np.log10(y[t + 1]) - np.log10(y[t])
    diff_x32 = np.log10(x[t + 1]) - np.log10(x[t])
    diff_y21 = np.log10(y[t]) - np.log10(y[t-1])
    diff_x21 = np.log10(x[t]) - np.log10(x[t-1])
    m = (.5*(diff_y32/diff_x32)) + (.5*(diff_y21/diff_x21))
    b = np.log10(y[t]) - m * np.log10(x[t])
    line = (10** b) * (x ** m)
    
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    axs.plot(x,y)
    axs.plot(x, line, label = "y = " + str(np.round(m,2)) + "x + " + str(np.round(b, 2)))
    axs.plot([x[t]], [y[t]],'o', 'k')
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.legend()
    axs.set_ylim(min(y), max(y))
    axs.set_aspect('equal')
    plt.show()
    
def cos(x, A, phi):
    period = 86400
    return A * np.cos(2 * np.pi * x / period  + phi)

def ResSquare(obs, exp):
    s = 0
    if(len(obs)!= len(exp)):
        print("The length of the arrays in the ResSquare function are not equal")
    for i in range(len(obs)):
        s += (obs[i] - exp[i])**2
    return s

def Variance(obs, exp, dof):
    return ResSquare(obs, exp)/dof

def Chi2(obs, exp, variance):
    return ResSquare(obs, exp)/variance

def Ftest(chi21, chi22, nu2):
    return nu2 * (chi21 - chi22) / chi22


def Harmonics(x, y, variance, GuessForAAndPhi = None): #Fitting function with F test included, this is only for Cosine  as of now.
    vars1, cov1 = curve_fit(cos, x, y, GuessForAAndPhi)
    period = 86400
    A = vars1[0]
    phi = vars1[1]
    w = 2 * np.pi /period
    def harm2(x, B):
        return(cos(x, A, phi) + B * np.cos(2*w*x + phi))
    vars2, cov2 = curve_fit(harm2, x, y)
    B = vars2[0]
    chi21 = Chi2(y, cos(x, A, phi), variance)
    chi22 = Chi2(y, harm2(x,B), variance)
    terms = 1
    nu2 = len(x) - terms - 1
    F = Ftest(chi21, chi22, nu2)
    print("F_1,2 = " + str(F) + "would you like to add another term? yes[y], no[n]")
    addTerm = raw_input()
    Y = y - harm2(x,B)
    Amplitudes = np.array([phi,A])
    zeros = np.zeros(len(Y))
    if(addTerm == 'y'):
        Amplitudes = np.append(Amplitudes, B)
        
    while addTerm == 'y':
        terms += 1
        def NewTerm(t, C):
            return C * np.cos((terms + 1) * w * t + phi)
        vars3, cov3 = curve_fit(NewTerm, x, Y)
        C = vars3[0]
        chi21 = Chi2(Y,zeros, variance)
        chi22 = Chi2(Y, NewTerm(x,C), variance)
        nu2 = len(x) - terms - 1
        F = Ftest(chi21, chi22, nu2)
        print("F_" + str(terms) + "," + str(terms + 1) + " = " + str(F) + "would you like to add another term? yes[y], no[n]")
        addTerm = raw_input()
        if (addTerm == 'y'):
            Amplitudes = np.append(Amplitudes, C)
        Y = Y - NewTerm(x, C)
    print("You entered something other than y, no more harmonics added. Last harmonic was number " + str(terms))
    print("All you need is to use the fuction BestFit and enter the coefficient array output by this function. The first coefficient is phi")
    return Amplitudes

def BestFit(x, Amplitudes): #Only use arrays output by the Harmonics function if you are unsure of how this function works. First element is phi, rest are amplitude coeff.
    Coeff = Amplitudes
    period = 86400
    w = 2 * np.pi / period
    phi = Amplitudes[0]
    Sum = Coeff[1] * np.cos(w * x + phi) 
    for i in range(len(Amplitudes) - 2):
        Sum += Coeff[i + 2] * np.cos((i + 2)*w*x + phi)
    return Sum

def Averages(moment_df, station0, binsize):
    """ Returns moment averages of neighbors of a certain station0
        Usually for stationary Trolley Runs
    Args:
        moment_df: ***An averaged moment dataframe (Moving average, using the pandas function)***, this works for multipole moments,
                   I am unsure asif it works for Legendre 2D moments. The format is the same as the
                   calc_moment_df (refer to analysis_helper) function that takes in an interpolation dataframe
        station0: The number of the "central station". Usually used for where the trolley
                  is underneath. 
        binsize: The size of the bins used in the moment_df (averaged moment df). The purpose of this is to delete the first trashy
                 entries of the arrays. Sorry, this was a poor way of doing this, I might fix this later.
    """
    #station0 should be an int
    nEvents = len(moment_df["st" + str(station0) + ",m1"].to_numpy()[(binsize-1):])
    avgsm1 = np.zeros((37, nEvents))
    avgsm2 = np.zeros((37, nEvents))
    avgsm3 = np.zeros((37, nEvents))
    avgsm4 = np.zeros((37, nEvents))
    avgsm5 = np.zeros((37, nEvents))
    
    avgms = np.zeros((5, 37, nEvents))
    
    avgsm1[0] = moment_df["st" + str(station0) + ",m1"].to_numpy()[(binsize-1):]
    avgsm2[0] = moment_df["st" + str(station0) + ",m2"].to_numpy()[(binsize-1):]
    avgsm3[0] = moment_df["st" + str(station0) + ",m3"].to_numpy()[(binsize-1):]
    avgsm4[0] = moment_df["st" + str(station0) + ",m4"].to_numpy()[(binsize-1):]
    avgsm5[0] = moment_df["st" + str(station0) + ",m5"].to_numpy()[(binsize-1):]
    
    for i in range(1,36):
        avgsm1[i] = (moment_df["st" + str((station0 + i)%72) + ",m1"].to_numpy()[(binsize-1):]\
        + moment_df["st" + str((station0 - i)%72) + ",m1"].to_numpy()[(binsize-1):])/2
        
        avgsm2[i] = (moment_df["st" + str((station0 + i)%72) + ",m2"].to_numpy()[(binsize-1):]\
        + moment_df["st" + str((station0 - i)%72) + ",m2"].to_numpy()[(binsize-1):])/2
        
        avgsm3[i] = (moment_df["st" + str((station0 + i)%72) + ",m3"].to_numpy()[(binsize-1):]\
        + moment_df["st" + str((station0 - i)%72) + ",m3"].to_numpy()[(binsize-1):])/2
        
        avgsm4[i] = (moment_df["st" + str((station0 + i)%72) + ",m4"].to_numpy()[(binsize-1):]\
        + moment_df["st" + str((station0 - i)%72) + ",m4"].to_numpy()[(binsize-1):])/2
        
        avgsm5[i] = (moment_df["st" + str((station0 + i)%72) + ",m5"].to_numpy()[(binsize-1):]\
        + moment_df["st" + str((station0 - i)%72) + ",m5"].to_numpy()[(binsize-1):])/2
        
        
    avgsm1[36] = moment_df["st" + str((station0 + 36)%72) + ",m1"].to_numpy()[(binsize-1):]
    
    avgsm2[36] = moment_df["st" + str((station0 + 36)%72) + ",m2"].to_numpy()[(binsize-1):]
    
    avgsm3[36] = moment_df["st" + str((station0 + 36)%72) + ",m3"].to_numpy()[(binsize-1):]
    
    avgsm4[36] = moment_df["st" + str((station0 + 36)%72) + ",m4"].to_numpy()[(binsize-1):]
    
    avgsm5[36] = moment_df["st" + str((station0 + 36)%72) + ",m5"].to_numpy()[(binsize-1):]
    
    return avgsm1,avgsm2,avgsm3,avgsm4,avgsm5

def SqRes(B, t, dt):
    """Square Residual method of quantifying noise. Better description to come
    
    This function is not used in current noise quantification as it is slow
    """
    N = int(len(B)/dt)
    s = 0
    for i in range(N):
        p = np.polyfit(t[i:i+dt], B[i:i+dt], 1)
        s+= np.sum(abs(B[i:i+dt] - np.polyval(p, t[i:i+dt]))**2)
    return s/len(t)

def FFTNoise(B, dt = 8):
    """FFT method of quantifying noise. Better description to come
    
    This function is not used in current noise quantification as is tricky to understand.
    """
    N = int(len(B)/dt)
    s = 0
    for i in range(N):
        C = np.fft.rfft(B[i:i+dt+1])
        s+= sum(abs(C[1:5]))
    return s