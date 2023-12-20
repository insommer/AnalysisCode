import matplotlib.pyplot as plt
import numpy as np 
#import matplotlib.patches as patches
#import lmfit
#from lmfit import Parameters
#import rawpy
#import imageio 
from scipy.optimize import curve_fit

kB = 1.38e-23
m = 9.988341e-27

def temperature_model(t, w0, T):
    
    #I define the constants explicitly since this function is passed to curve fit
    kB = 1.38e-23 #Boltzmann's constant
    m = 9.988341e-27 #Li-6 mass in kg    
    t0 = 0
    
    model = w0*np.sqrt(   1 +       (kB/m)*abs(T)*(t-t0)**2/(w0**2)   )
    # model = w0*np.sqrt((kb*T*(t-t0)**2)/(m*w0**2))
    return model



def temperature_fit(params, widths_array, tof_array):
    
    #fit_array = np.linspace(tof_array[0], tof_array[-1], len(tof_array))
    min_time = min(tof_array)
    max_time = max(tof_array)
    
    min_width = min(widths_array)
    max_width = max(widths_array)
    
    # initial guess for parameters of temperature_model
    w0guess = min_width
    slope = (max_width-min_width)/(max_time-min_time)
    ##Tguess = (slope)**2*params.m/params.kB
    Tguess = (slope)**2 * m/kB
    
    # returns optimal values for fitting popt, and covariance array pcov
    popt, pcov = curve_fit(temperature_model, tof_array, widths_array, p0 = [w0guess, Tguess])
    
    # arrays to plot fitted function
    times_fit = np.linspace(min_time, max_time, 100)
    widths_fit = temperature_model(times_fit, popt[0], popt[1])
    
    return tof_array, times_fit, widths_fit, popt, pcov

 
def thermometry(params, images, tof_array, do_plot = False, data_folder = None):
    
    widthsx = np.zeros(len(images))
    widthsy = np.zeros(len(images))
    
    #fill arrays for widths in x and y directions
    for index, image in enumerate(images):
        widthsx[index], x, widthsy[index], y  = fitgaussian(image)
        widthsx[index] = widthsx[index]*params.camera.pixelsize_meters/params.magnification
        widthsy[index] = widthsy[index]*params.camera.pixelsize_meters/params.magnification
        if index == 0:
            print("widthx = "+str(widthsx[index]*1e6)+" um")
            print("widthy = "+str(widthsy[index]*1e6)+" um")
    #these plots will still show even if the fit fails, but the plot underneath the fit will not     
    # if (do_plot):
    #     plt.figure()
    #     plt.xlabel("Time of flight (ms)")
    #     plt.ylabel("1/e^2 width x of atom cloud (uncalibrated units)")
    #     plt.scatter(tof_array, widthsx)
        
    #     plt.figure()
    #     plt.xlabel("Time of flight (ms)")
    #     plt.ylabel("1/e^2 width y of atom cloud (uncalibrated units)")
    #     plt.scatter(tof_array, widthsy)        
        
    fitx_array, plotx_array, fitx, poptx, pcovx = temperature_fit(params, widthsx, tof_array)
    fity_array, ploty_array, fity, popty, pcovy = temperature_fit(params, widthsy, tof_array)
    if (do_plot):
        #plot the widths vs. position along x direction
        plt.figure()
        plt.title("Temperature fit x, T = {} uK".format(poptx[1]*1e6))
        plt.xlabel("Time of flight (ms)")
        plt.ylabel("width of atom cloud (um)")
        plt.scatter(1e3*tof_array, 1e6*widthsx)
        plt.plot(1e3*plotx_array, 1e6*temperature_model(plotx_array, *poptx))   
        if data_folder:
            plt.savefig(data_folder+r'\\'+"temperature x.png", dpi = 500)
        #plot the widths vs. position along y direction
        plt.figure()
        plt.title("Temperature Fit y, T = {} $\mu$K".format(popty[1]*1e6))
        plt.xlabel("Time of Flight (ms)")
        plt.ylabel("Width of Atom Cloud ($\mu$m)")
        plt.scatter(1e3*tof_array, 1e6*widthsy)
        plt.plot(1e3*ploty_array, 1e6*temperature_model(ploty_array, *popty)) 
        if data_folder:
            plt.savefig(data_folder+r'\\'+"temperature y.png", dpi = 500)        
    return poptx, pcovx, popty, pcovy
   
def thermometry1D(params, columnDensities, tof_array, thermometry_axis="x", 
                  do_plot = False, save_folder = None, reject_negative_width=False,
                  newfig=True):
    #1. Find cloud size (std dev) vs time
    widths=[]
    times=[]
    numbers=[]
    dx = params.camera.pixelsize_meters/params.magnification
    for index, density2D in enumerate(columnDensities):
        
        density1D = integrate1D(density2D, dx=dx, free_axis=thermometry_axis)
        
        xdata = np.arange(np.shape(density1D)[0])*dx
        
        popt_gauss  = fitgaussian1D(density1D,xdata,dx=dx, doplot=True, xlabel=thermometry_axis, ylabel="density")
        
        
        if popt_gauss is not None: #fit succeeded
            if popt_gauss[2] >0 or not reject_negative_width:
                w = abs(popt_gauss[2])
                widths.append(w)
                times.append(tof_array[index])
                numbers.append(popt_gauss[0]* w*(2*np.pi)**0.5 ) #amp  = N/(w*(2*np.pi)**0.5)
    
    numbers=np.array(numbers)
    
    widths=np.array(widths)
    
    times=np.array(times)       
    #2. Fit to a model to find temperature    
    # fitx_array, plotx_array, fit, popt, pcov = temperature_fit(params, widths, times)
    try:
        tof_array, times_fit, widths_fit, popt, pcov = temperature_fit(params, widths, times)
    except RuntimeError:
        popt = None
        pcov = None
        
    if (do_plot):
        #plot the widths vs. time
        if newfig:
            plt.figure()
        plt.rcParams.update({'font.size': 14})
        AxesAndTitleFont = 20
        if popt is not None:
            plt.title("{0}: T = {1:.2f} $\mu$K".format(thermometry_axis, popt[1]*1e6), fontsize = AxesAndTitleFont)
            plt.plot(1e3*times_fit, 1e6*widths_fit, color = 'blue', zorder =1)
        plt.xlabel("Time of Flight (ms)", fontsize = AxesAndTitleFont)
        plt.ylabel("Std. dev ($\mu$m)", fontsize = AxesAndTitleFont)
        plt.scatter(tof_array/1e-3, widths/1e-6, color = 'red', zorder = 2)
        
        
        plt.tight_layout()
        if save_folder:
            plt.savefig(save_folder+r'\\'+"temperature {}.png".format(thermometry_axis), dpi = 300)
        # plt.figure()
        # plt.plot(1e3*tof_array, numbers,'o')
        # plt.xlabel("Time of flight (ms)")
        # plt.ylabel("Atom number")
        # plt.title("Atom Number {}".format(thermometry_axis))
        plt.tight_layout()
        if save_folder:
            plt.savefig(save_folder+r'\\'+"atom number {}.png".format(thermometry_axis), dpi = 300)
    return popt, pcov


def Gaussian(x, amp, center, w, offset):
    #amp = N/(w*(2*np.pi)**0.5)
    return amp*np.exp(-.5*(x-center)**2/w**2) + offset


   
def exponential(x, a, tau, c):
    return a * np.exp(-x/tau) + c    


    
def integrate1D(array2D, dx=1, free_axis="y"):
    #free_axis is the axis that remains after integrating
    if free_axis == 'x':
        axis = 0
    elif free_axis == 'y':
        axis = 1
    array1D = np.sum(array2D, axis = axis)*dx
    return array1D


def fit_exponential(xdata, ydata ,dx=1, doplot = False, label="", title="", 
                    newfig=True, xlabel="",ylabel="", offset = None, legend=False):
    

    #fit for the parameters a , b, c
    a = max(ydata) - min(ydata)
    tau = (max(xdata)-min(xdata))/2
    c = min(ydata)
    xfit = np.linspace(min(xdata),max(xdata), 1000)
    
    if offset is None:
        func = exponential
        guess= [a,tau,c]
        label = 'fit: a=%5.3f, tau=%5.3f, c=%5.3f'
    else:
        func = lambda x,a,tau: exponential(x,a,tau,offset)
        guess = [a,tau]
        label = 'fit: a=%5.2e\n tau=%5.2e\n c={:.2e} (Fixed)'.format(offset)
        
    popt, pcov = curve_fit(func, xdata, ydata, p0=guess)       

    #poptarray([2.56274217, 1.37268521, 0.47427475])
    plt.plot(xfit, func(xfit, *popt), 'r-', label= label % tuple(popt))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend:
        plt.legend()
    plt.show()
    return popt, pcov

def fitgaussian(array, do_plot = False, vmax = None,title="", 
                save_column_density = False, column_density_xylim = None): 
    #np.sum(array, axis = 0) sums over rows, axis = 1 sums over columns
    rows = np.linspace(0, len(array), len(array))
    cols = np.linspace(0, len(array[0]), len(array[0]))
    row_sum = np.sum(array, axis = 0)  
    col_sum = np.sum(array, axis = 1)
    # print("rows = "+str(rows))
    # print("np.shape(array) = "+str(np.shape(array)))
    # print("np.shape(rows) = "+str(np.shape(rows)))
    # print("np.shape(cols) = "+str(np.shape(cols)))
    # print("np.shape(row_sum) = "+str(np.shape(row_sum)))
    # print("np.shape(col_sum) = "+str(np.shape(col_sum)))
    ampx = np.max(row_sum)
    centerx = np.argmax(row_sum)
    wx = len(rows)/12
    # offsetx = row_sum[0]
    ampy = np.max(col_sum)
    centery = np.argmax(col_sum)
    wy = len(cols)/12
    # offsety = col_sum[0]
    # if (do_plot):
    #     plt.figure()
    #     plt.imshow(array, cmap = 'gray')
    #     plt.figure()
    #     plt.title("col sum (for x direction fit)")
    #     plt.plot(cols, row_sum)
    #     plt.xlabel("pixel index")
    #     plt.ylabel("sum over array values")
    #     plt.figure()
    #     plt.title("row sum (for y direction fit)")
    #     plt.plot(rows, col_sum)
    #     plt.xlabel("pixel index")
    #     plt.ylabel("sum over array values")  
    
    if do_plot:
        #see the input array
        plt.rcParams.update({'font.size' : 10})
        plt.figure(figsize=(12,5))
        plt.subplot(121)
        if vmax == None:
            vmax = array.max()
        
        plt.imshow(array, cmap = 'jet',vmin=0,vmax=vmax)
        
        if column_density_xylim == None:
            column_density_xylim = np.zeros(4)
            # column_density_xylim[1] = len(array[0])
            # column_density_xylim[3] = len(array)
            column_density_xylim[1], column_density_xylim[3] = array.shape[1::-1]
            
        print("column_density_xylim = "+str(column_density_xylim))
        column_density_xylim = np.array(column_density_xylim)
        if column_density_xylim[1] == -1:
                column_density_xylim[1] = len(array[0])
        if column_density_xylim[3] == -1:
                column_density_xylim[3] = len(array) 
        plt.xlim(column_density_xylim[0], column_density_xylim[1])
        plt.ylim(column_density_xylim[3], column_density_xylim[2])
        plt.title(title)
        plt.colorbar(pad = .1)
        if save_column_density:
            plt.savefig(title + ".png", dpi = 600)
        #plot the sum over columns
        #plt.figure()
        plt.subplot(122)
        #plt.title("col sum (fit in x direction)")
        plt.plot(cols, row_sum, label="data_vs_x")
        
        plt.xlabel("pixel index")
        plt.ylabel("sum over array values")
        #plot the sum over rows
        #plt.figure()
        #plt.title("row sum (fit in y direction)")
        plt.plot(rows, col_sum, label="data vs y")
        
        plt.xlabel("pixel index")
        plt.ylabel("sum over array values")
        plt.tight_layout()
        plt.legend()
    
        plt.plot(cols, Gaussian(cols, *[ampx, centerx, wx,0]), label="guess vs x")
        plt.plot(rows, Gaussian(rows, *[ampy, centery, wy,0]), label="guess vs y")  
        plt.legend()
        plt.tight_layout()
        
    widthx, center_x, widthy, center_y = np.nan, np.nan, np.nan, np.nan
    try:
        poptx, pcovx = curve_fit(Gaussian, cols, row_sum, p0=[ampx, centerx, wx,0])
        widthx = abs(poptx[2])
        center_x = poptx[1]
        
        if (do_plot):
            plt.plot(cols, Gaussian(cols, *poptx), label="fit vs x")
            plt.legend()
            plt.tight_layout()
    
    except RuntimeError as e:
        print(e)
        
    try:
        popty, pcovy = curve_fit(Gaussian, rows, col_sum, p0=[ampy, centery, wy,-1e13])
        widthy = abs(popty[2])
        center_y = popty[1]  
        
        if (do_plot):
            plt.plot(rows, Gaussian(rows, *popty), label="fit vs y")  
            plt.legend()
            plt.tight_layout()
    
    except RuntimeError as e:
        print(e)
        
    return widthx, center_x, widthy, center_y


def fitgaussian1D(data , xdata=None, dx=1, doplot = False, 
                  subtract_bg=True, signal_feature='narrow', 
                  label="", title="", newfig=True, 
                  xlabel="", ylabel="", xscale_factor=1, legend=False,
                  yscale_factor=1):
    
    if subtract_bg:
        bg = fitbg(data, signal_feature=signal_feature) 
        originalData = data.copy()
        data = data - bg
        
        offset_g = 0
    else:
        offset_g = offset_g = min( data[:10].mean(), data[-10:].mean() )
    
    datalength = len(data)
    
    if xdata is None:
        xdata = np.arange( datalength )*dx  
        
    #initial guess:
    amp_g = data.max()
    center_g = xdata[ data.argmax() ]    
    w_g = ( data > 0.6*data.max() ).sum() * dx / 2
    
    guess = [amp_g, center_g, w_g, offset_g]
    
    try:
        popt, pcov = curve_fit(Gaussian, xdata, data, p0=guess, bounds=([-np.inf, -np.inf, 0, -np.inf],[np.inf]*4) )
    except Exception as e:
        print(e)
        return None  
 
    #      
    if doplot:
        if newfig:
            plt.figure()            
        if subtract_bg:                
            plt.plot(xdata*xscale_factor, originalData*yscale_factor, '.', label="{} data".format(label))
            plt.plot(xdata*xscale_factor, (Gaussian(xdata,*popt)+bg) * yscale_factor, label="{} fit".format(label))
            plt.plot(xdata*xscale_factor, bg*yscale_factor, '.', markersize=0.3)
            # ax.plot(xdata*xscale_factor, (Gaussian(xdata,*guess)+bg) * yscale_factor, label="{} fit".format(label))
        else:
            plt.plot(xdata*xscale_factor, data*yscale_factor, '.', label="{} data".format(label))
            plt.plot(xdata*xscale_factor, Gaussian(xdata,*popt) * yscale_factor, label="{} fit".format(label))

    if doplot:
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if legend:
            plt.legend()
    return popt




def fitbg(data, signal_feature='narrow', signal_width=10, fitbgDeg=5): 
       
    datalength = len(data)
    signalcenter = data.argmax()
    datacenter = int(datalength/2)    
    xdata = np.arange(datalength)
    
    if signal_feature == 'wide':
        mask_hw = int(datalength/3)
        bg_mask = np.full(xdata.shape, True)
        bg_mask[signalcenter - mask_hw: signalcenter + mask_hw] = False  
        
        p = np.polyfit( xdata[bg_mask], data[bg_mask], deg=2 )        
        
    else:
        mask_hw = int(datalength/signal_width)
        bg_mask = np.full(xdata.shape, True)
        center_mask = bg_mask.copy()
        bg_mask[signalcenter - mask_hw: signalcenter + mask_hw] = False  
        center_mask[datacenter - mask_hw : datacenter + mask_hw] = False        
        bg_mask = bg_mask * center_mask
        bg_mask[:mask_hw] = True
        bg_mask[-mask_hw:] = True
        
        p = np.polyfit( xdata[bg_mask], data[bg_mask], deg=fitbgDeg )
    
    return np.polyval(p, xdata)
