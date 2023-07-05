from ImageAnalysis import ImageAnalysisCode
from ImageAnalysis import SIUnits as units
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from scipy.optimize import curve_fit


def AnalyzeZyla(data_folder, rowstart = 0, rowend =-1, columnstart = 0, 
                columnend = -1,picturesPerIteration = 3, iterationNum=None, verbose=True, angle_deg=1):
    
    params = ImageAnalysisCode.ExperimentParams(picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
    if iterationNum is None:
        images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)
    else:
        images_array = ImageAnalysisCode.LoadFromSpooledSeries(params, iterationNum, data_folder=data_folder)
    
    # ImageAnalysisCode.ShowImagesTranspose(images_array)
    
    Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                    firstFrame=0, correctionFactorInput=1.0,  
                    subtract_burntin=0, preventNAN_and_INF=False)
    # plt.figure()
    # plt.imshow(np.array(images_array[0][0]-images_array[0][2],dtype=np.float64)/(images_array[0][1]-images_array[0][2]),vmin=0,vmax=1.1)
    # plt.imshow(images_array[0][0]-images_array[0][1])
    
    imgNo = len(columnDensities)
    #angle_deg= 1 #rotates ccw
    
    AtomNumberList=[]
    widths_x = []
    widths_y = []
    
    for ind in range(imgNo):
        rotated_ = rotate(columnDensities[ind], angle_deg, reshape = False)[rowstart:rowend,columnstart:columnend]
        # rotated_=columnDensities[ind]
        if ind==0: #first time
            rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
        rotated_columnDensities[ind] = rotated_
    
        #preview:
        dx=params.camera.pixelsize_meters/params.magnification
        popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = True, title="1D density",
                                                      ylabel1D="1d density (atoms/$\mu$m)", xlabel1D="position ($\mu$m)",
                                                      title2D="column density",
                                                      xscale_factor=1/units.um, yscale_factor=units.um)
        
        wy = abs(popt1[2])
        AtomNumberY = popt1[0]* wy*(2*np.pi)**0.5 
        AtomNumberList.append(AtomNumberY)
        print("\n{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
        # print(popt1)
        
        if popt0 is not None:
            print("RMS cloud size x: {:.2f} um".format(popt0[2]/units.um))
            print("RMS cloud size y: {:.2f} um".format(popt1[2]/units.um))
        
            widths_x.append(popt0[2])
            widths_y.append(popt1[2])
    
    if verbose:
        print('\nThe average number of atoms:{:.2e}'.format(np.mean(AtomNumberList)))        
        print("Mean RMS width x: {:.2f} +/- {:.2f} um".format(np.mean(widths_x)/units.um, np.std(widths_x)/units.um))
        print("Mean RMS width y: {:.2f} +/- {:.2f} um".format(np.mean(widths_y)/units.um, np.std(widths_y)/units.um))
    
def AnalyzeChameleonBinary(data_folder, rowstart = 0, rowend =-1, columnstart = 0, 
                columnend = -1,picturesPerIteration = 3, t_exp = 10e-6, binsize = 1, iterationNum=None, CircularMask = None, centerx = 626, centery = 580,
                radius = 260):
    
           params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "chameleon")      
           # images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)
           # images_array = ImageAnalysisCode.loadSeriesRAW(params = params, picturesPerIteration=picturesPerIteration, data_folder = data_folder)
           images_array = ImageAnalysisCode.loadSeriesPGM(picturesPerIteration=picturesPerIteration, data_folder = data_folder, 
                                                          binsize=binsize, file_encoding = 'binary')



           # Number_of_atoms, columnDensities = ImageAnalysisCode.flsImaging(images_array, params = params, firstFrame=0, rowstart = 0, rowend = -1, 
           #                                                                 columnstart =0, columnend = -1, subtract_burntin = 1)


           ImageAnalysisCode.ShowImagesTranspose(images_array)


           # print(data_folder+r'\\image-02242023185744-1.Raw')
           # print(ImageAnalysisCode.loadRAW(data_folder+'//image-02242023185744-1.Raw'))

           Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                           firstFrame=0, correctionFactorInput=1, rowstart = rowstart, rowend = rowend, columnstart = columnstart, columnend = columnend, 
                           subtract_burntin=0, preventNAN_and_INF=True)
           if CircularMask:
               columnDensities = ImageAnalysisCode.CircularMask(columnDensities, centerx=centerx/binsize, centery=centery/binsize, 
                                                                radius=radius/binsize)
           plt.figure()
           plt.imshow(ratio_array[0], vmin = 0, vmax=1.5, cmap = 'gray')
           plt.colorbar()
           plt.show()


           # for count, x in enumerate(columnDensities):       
           #     widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count], do_plot = 1)
           #     print("Center x:",center_x)
           #     print("Center y:",center_y)
           # tof_array = np.loadtxt(data_folder + '/tof_list.txt')*1e-3+t0
           # poptx, pcovx, popty, pcovy = ImageAnalysisCode.thermometry(params, columnDensities, tof_array, do_plot = True, data_folder = data_folder)


           center_x_array = np.zeros(len(images_array))
           center_y_array = np.zeros(len(images_array))
           for count, x in enumerate(columnDensities):       
               widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count],title = "Vertical Column Density",\
                                                                                   vmax = None, do_plot = 1, save_column_density=0,\
                                                                                       column_density_xylim=(columnstart, columnend, rowstart, rowend))
               center_x_array[count] = center_x
               center_y_array[count] = center_y
               print("Center x:",center_x)
               print("Center y:",center_y)
               print("Number of atoms:{}e6".format(round(Number_of_atoms[count]/(1e6))))
               if (count+1)%2 == 0:
                   print("difference in center x:",center_x_array[count]-center_x_array[count - 1])
                   print("difference in center y:",center_y_array[count]-center_y_array[count - 1])

           '''
           ImageAnalysisCode.imageFreqOptimization(np.loadtxt(data_folder+"/imgfreq.txt"), Number_of_atoms, ratio_array)
           plt.imshow(ratio_array[0][rowstart:rowend,columnstart:columnend],vmin=0,vmax=1.2,cmap="gray")
           densityvsrow = np.sum(n2d[0][rowstart:rowend,columnstart:columnend], 1)
           print("densityvsrow = "+str(np.shape(densityvsrow)))
           plt.figure(figsize=(4,3))
           plt.plot(densityvsrow)
           '''

           def gaussianBeam(x, amp, center, w, offset, slope):
               return offset + amp*np.exp(-2*(x-center)**2/w**2) + slope*x


           def fitgaussian(xdata, ydata, do_plot=True):
               popt, pcov = curve_fit(gaussianBeam, xdata, ydata,p0=[3e9, 670, 10, 5e9, 3e7])
               #print(popt)
               if (do_plot):
                   plt.plot(xdata,ydata)
                   plt.plot(xdata,gaussianBeam(xdata, *popt))
                   plt.title("amplitude = {:.2e}".format(popt[0]))
                   plt.ylabel("1D atomic density arb units")
                   plt.xlabel("vertical row index")
                   plt.tight_layout()
                   # plt.savefig("atom count plot.png", dpi = 600)
               return popt,pcov

           '''
           #odt-specific code below
           angle_deg= 0 #rotates ccw
           rotated_columnDensities = rotate(columnDensities[0][rowstart:rowend,columnstart:columnend], angle_deg, reshape = False)
           plt.imshow(rotated_columnDensities)
           plt.colorbar()
           densityvsrow = np.sum(rotated_columnDensities, 1)*deltaY
           print("densityvsrow = "+str(np.shape(densityvsrow)))
           plt.figure(figsize=(4,3))
           plt.plot(densityvsrow)

           nstart = 659
           nstop = 685
           xvalues = np.arange(nstart, nstop)
           popt, pcov = fitgaussian(xvalues, densityvsrow[nstart:nstop], do_plot = 1)
           odt_fit = gaussianBeam(xvalues, popt[0], popt[1], popt[2], popt[3], popt[4])-popt[4]*xvalues
           num_atoms = simpson(odt_fit, xvalues)*deltaX
           print("Number of atoms in ODT: {}e6".format(num_atoms/(1e6)))
           plt.show()  
           '''
        