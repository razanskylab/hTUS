#-*- coding: utf-8 -*-
# Some basic tools for data processing
from numpy import nonzero, sign
from scipy import stats
from scipy.interpolate import interp1d
from pylab import *
from scipy import signal
from scipy.fft import fft, ifft, fftshift
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors
import colormaps as cmaps
from math import factorial
import getpass
 
 
def new_rgba_colormap(vrgba,Np,outFlag):
  # Builds colormap by means of linear interpolation between the points of the
  # Nx5 array vrgba (value,red, green, blue, opacity), where value lies between
  # 0 and 1.
  # Np is the total number of points to perform the interpolation.
  # If outFlag=0 it returns a colormap object, else it returns a rgba array.
  # It distributes the colors evenly in the range (0,1).

  lims=squeeze(vrgba[:,0])
  redy=squeeze(vrgba[:,1])
  greeny=squeeze(vrgba[:,2])
  bluey=squeeze(vrgba[:,3])
  opcty=squeeze(vrgba[:,4])

  nLevels=len(lims)
  x=linspace(0,1,Np)
  green=zeros(Np)
  red=zeros(Np)
  blue=zeros(Np)
  alpha=zeros(Np)
  gForce=zeros((Np,4))

  for nn in range(nLevels-1):
    inSide=nonzero(logical_and(lims[nn]<=x, x <=lims[nn+1]))[0]
    if inSide.any() :
      red[inSide]=(redy[nn+1]-redy[nn])/(lims[nn+1]-lims[nn])*(x[inSide]-lims[nn])+redy[nn]
      green[inSide]=(greeny[nn+1]-greeny[nn])/(lims[nn+1]-lims[nn])*(x[inSide]-lims[nn])+greeny[nn]
      blue[inSide]=(bluey[nn+1]-bluey[nn])/(lims[nn+1]-lims[nn])*(x[inSide]-lims[nn])+bluey[nn]
      alpha[inSide]=(opcty[nn+1]-opcty[nn])/(lims[nn+1]-lims[nn])*(x[inSide]-lims[nn])+opcty[nn]
    else:
      print('Not enough points to define the colormap')
      return

  gForce[:,0]=red; gForce[:,1]=green; gForce[:,2]=blue; gForce[:,3] = alpha;
  kk=nonzero(gForce>1)[0];
  gForce[kk]=1.0;
  kk=nonzero(gForce < 0)[0];
  gForce[kk]=0.0;
  ## Debugging  
  #figure()
  #plot(gForce)
  #show()
  #print(repr(shape(gForce))+' '+repr(shape(alpha)))
  if (outFlag):
    newCmap=gForce
  else:
    cmapName='rgbaColormap'
    newCmap = matplotlib.colors.ListedColormap(gForce, cmapName, N=None)
    newCmap._init()
    #print(repr(shape(newCmap._lut)))
    #newCmap._lut[0:len(alpha),3]=alpha
    #figure()
    #plot(newCmap._lut)
    #show()

  return newCmap

def rgbaCmap(cmapName,alpha,outFlag):
  # Given the name of a colormap and the alpha
  # channel vector, it returns an rgba
  # version of the colormap. If outFlag=0 it returns
  # a colormap object, else it returns a rgba array.
  cdf=linspace(0.0,1.0,len(alpha))
  gForce=eval('cm.'+cmapName+'(cdf)')
  #print(repr(shape(gForce))+' '+repr(shape(alpha)))
  gForce[:,3] = alpha
  if (outFlag):
    newCmap=gForce
  else:
    cmapName=cmapName+'Trans'
    newCmap = matplotlib.colors.ListedColormap(gForce,cmapName,N=None)
    newCmap._init()
    #print(repr(shape(newCmap._lut)))
    newCmap._lut[0:len(alpha),3]=alpha
    #figure()
    #plot(newCmap._lut)
    #show()

  return newCmap

def pageSetup(fig_width_mm,aratio):
  # Fine tuning for printing stuff ======================================
  fig_width_pt = 2.834645669*fig_width_mm   # Get this from LaTeX using \showthe\columnwidth
  golden_mean = aratio#1.0/sqrt(2) #(sqrt(5)-1.0)/2.0         # Aesthetic ratio

  inches_per_pt = 1.0/72.27               # Convert pt to inch
  fig_width = fig_width_pt*inches_per_pt  # width in inches
  fig_height = fig_width*golden_mean      # height in inches
  fig_size =  [fig_width,fig_height]

  if getpass.getuser() == "ektoras":
      params = {'backend': 'pdf',
                'axes.labelsize': 10,
                'axes.linewidth': 0.8,
                'font.size': 10,
                'font.family': 'sans-serif',
                # 'text.font.family': 'sans-serif',
                # 'axes.font.family': 'sans-serif',
                'boxplot.flierprops.linewidth': 1.0,
                'boxplot.boxprops.linewidth' : 1.0,
                'legend.fontsize': 'small',
                'xtick.labelsize': 8,
                'ytick.labelsize': 8,
                'text.usetex': False,
                'figure.figsize': fig_size}

  if getpass.getuser() == "hofmannu":
     params = {'backend': 'pdf',
               'axes.labelsize': 10,
               'axes.linewidth': 0.8,
               'font.size': 10,
               'font.family': 'sans-serif',
               'boxplot.flierprops.linewidth': 1.0,
               'boxplot.boxprops.linewidth' : 1.0,
               'legend.fontsize': 'small',
               'xtick.labelsize': 8,
               'ytick.labelsize': 8,
               'figure.figsize': fig_size}
  else:
    params = {'backend': 'pdf',
                'axes.labelsize': 10,
                'axes.linewidth': 0.8,
                'font.size': 10,
                'font.family': 'sans-serif',
                # 'text.font.family': 'sans-serif',
                # 'axes.font.family': 'sans-serif',
                'boxplot.flierprops.linewidth': 1.0,
                'boxplot.boxprops.linewidth' : 1.0,
                'legend.fontsize': 'small',
                'xtick.labelsize': 8,
                'ytick.labelsize': 8,
                'text.usetex': False,
                'figure.figsize': fig_size}    
   
  rcParams.update(params)
  return


def tukeywin(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.

    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output

    Reference
    ---------

http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html

    '''
    # Special cases
    if alpha <= 0:
        return ones(window_length) #rectangular window
    elif alpha >= 1:
        return hanning(window_length)

    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)

    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + cos(2*pi/alpha * (x[first_condition] - alpha/2) ))

    # second condition already taken care of

    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + cos(2*pi/alpha * (x[third_condition] - 1 + alpha/2)))

    return w

