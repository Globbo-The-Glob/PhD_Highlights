import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 

def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])
  
def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def graph_tendon(radius,S,L,steps,file):
  #No. of blocks
  ## Splines
  shores = np.linspace(15,60,steps)
  colour = linear_gradient("#83f5e5", finish_hex="#e761bd", n=steps)
  indexs = np.zeros((1,steps))
  for i in range(1,steps):
    indexs[0,i] = find_nearest(shores, S[i]) # pulls index for colour
  
  fig, ax = plt.subplots()
  coloursrgb = np.zeros((steps,3))
  coloursrgb[:,0] = colour['r'][:]
  coloursrgb[:,1] = colour['g'][:]
  coloursrgb[:,2] = colour['b'][:]
  cmap_name = 'my_list'
  map = LinearSegmentedColormap.from_list(cmap_name, coloursrgb/255, N=steps)
      
  for step in range(1,steps):
      # Draw polygon 
      # Edges of polygon in general are spline points 
      coord = [[(L/steps)*step,radius[step]], [(L/steps)*(step),-radius[step]], [(L/steps)*(step-1),-radius[step-1]], [(L/steps)*(step-1),radius[step-1]]]
      coord.append(coord[0]) #repeat the first point to create a 'closed loop'
      xs, ys = zip(*coord) #create lists of x and y values
      ii = int(indexs[0,step])
      ax.fill(ys,xs,colour['hex'][ii])
      ax.set_xlabel('x [m]')
      ax.set_ylabel('y [m]')
      ax.set_xlim([-0.05, 0.05])    
  fig.suptitle(file)
  fig.savefig(file + '.png')

def compare_coarse(params,steps,file):
  
  s1 = params[0]
  s2 = params[1]
  s3 = params[2]
  s4 = params[3]
  r1 = params[4]
  r2 = params[5]
  L = params[6]
  
  x = np.linspace(0,1,steps) 
  S = s1*(1-x)**3 + s2*3*x*(1-x)**2 + s3*3*x**2*(1-x) + s4*x**3
  radius = 0.01*(1-x)**3 + r1*3*x*(1-x)**2 + r2*3*x**2*(1-x) + 0.01*x**3
  graph_tendon(radius,S,L,steps,file+"coarse")
  
  x = np.linspace(0,1,1000) 
  S = s1*(1-x)**3 + s2*3*x*(1-x)**2 + s3*3*x**2*(1-x) + s4*x**3
  radius = 0.01*(1-x)**3 + r1*3*x*(1-x)**2 + r2*3*x**2*(1-x) + 0.01*x**3
  graph_tendon(radius,S,L,1000,file+"fine")
