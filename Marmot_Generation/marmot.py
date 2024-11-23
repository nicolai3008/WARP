import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
import warnings
from tqdm import tqdm
import pandas as pd
from scipy.optimize import minimize
import datetime
import circle_fit as cf
import shutil
from subprocess import call
warnings.filterwarnings("ignore")
from mcstas_functions import *
import math

# Base functions
def E2lambda(E):
    return 9.045/np.sqrt(E)

def lambda2E(lam):
    return (9.045/lam)**2

def mcstas_x_E_plot(folder):
    file = "E_end.x_E"
    intensity, error, count, parameters = read_data_files(folder, file, D=2)
    plot_heatmap(intensity, parameters, save_folder=folder, save_name="x_E")
    
def mcstas_resolution_function(folder, d_pars, E_pars):
    file = "E_end.x_E"
    intensity, error, count, parameters = read_data_files(folder, file, D=2)
    x, E = xy(parameters,intensity)
    E = E[:,0]
    averages = []
    stds = []
    for i in range(len(x[0])):
        i_row = intensity[:,i]
        try:
            average = np.average(E,weights=i_row)
            std = resolution(i_row,E)*2.355
            averages.append(average)
            stds.append(std)
        except ZeroDivisionError:
            averages.append(np.nan)
            stds.append(np.nan)
        
    averages = np.array(averages)
    stds = np.array(stds)
    
    # Find a way to flip if necessary
    
    # Ideal resolution function
    xs = np.linspace(min(x[0]),max(x[0]),1000)
    xss = np.linspace(0,1,1000)
    Es = energy_det(xss,E_pars)
    xs = xs[::-1]
    
    plt.figure()
    plt.plot(x[0],averages)
    plt.plot(xs,Es,'k--')
    plt.fill_between(x[0], averages-stds, averages+stds, alpha=0.5)
    plt.xlabel('x / m')
    plt.ylabel('Energy / meV')
    plt.title('Resolution function')
    plt.savefig(f'{folder}/resolution.png')
    plt.close()
    
    plt.figure()
    plt.plot(averages,stds)
    plt.xlabel('Energy / meV')
    plt.ylabel('Resolution / meV')
    plt.title('Resolution function')
    plt.savefig(f'{folder}/resolution2.png')
    plt.close()

def circle(x,y,R):
    theta = np.linspace(0,2*np.pi,10000)
    x_c = x + R*np.cos(theta)
    y_c = y + R*np.sin(theta)
    return x_c, y_c

def circle_in_area(x,y,R,angle=2):
    x_c, y_c = circle(x,y,R)
    idxs = np.logical_and(x_c < 5, x_c > 0.5)
    x_c = x_c[idxs]
    y_c = y_c[idxs]
    idxs = np.logical_and(y_c < np.tan(angle*np.pi/360)*x_c,y_c > -np.tan(angle*np.pi/360)*x_c)
    x_c = x_c[idxs]
    y_c = y_c[idxs]
    return x_c, y_c

def x_0(R_x,R_y,R):
    x = R_x + np.sqrt(R**2 - R_y**2)
    return x


# Dist -> Angles in space -> Marmot curves
def gen_folder(folder):
    with open(f"{folder}/parameter.txt", "r") as f:
        lines = f.readlines()
        if len(lines) == 2:
            n = 0
        else:
            for line in lines[::-1]:
                line = line.split(";")
                if line[0] == "\n":
                    continue
                n = int(line[0])
                if n != 0:
                    break
            
    n += 1
    
    
    return n, f'{folder}/Data/{n}'
    
def points_in_space(mm=0.001, angle=4):
    Nx = int((5-0.25)//mm)
    Ny = int(2*np.tan(angle*np.pi/360)*5//mm)
    x_space = np.linspace(0.25,5,Nx)
    y_space = np.linspace(-np.tan(angle*np.pi/360)*5,np.tan(angle*np.pi/360)*5,Ny)
    
    return x_space, y_space

def energy_det(x_per,par):
    # Fix to work with any distribution
    #x_per = 1-x_per
    E_0, E_1, dE_0, dE_1, ddE_0 = par
    a0 = E_0
    a1 = dE_0
    a2 = ddE_0/2
    a3 = -3*dE_0-4*E_0+4*E_1-ddE_0-dE_1
    a4 = 2*dE_0+3*E_0-3*E_1+ddE_0/2+dE_1
    E = a0 + a1*x_per + a2*x_per**2 + a3*x_per**3 + a4*x_per**4
    return E

def inv_energy_det(E,par):
    E_0, E_1, dE_0, dE_1, ddE_0 = par
    a0 = E_0
    a1 = dE_0
    a2 = ddE_0/2
    a3 = -3*dE_0-4*E_0+4*E_1-ddE_0-dE_1
    a4 = 2*dE_0+3*E_0-3*E_1+ddE_0/2+dE_1
    x = np.roots([a4,-a3,a2,-a1,a0-E])
    x = x[np.isreal(x)]
    x = x[np.logical_and(x>=0,x<=1)]
    return x[0]

def pos_det(x_per,detector):
    x = x_per*detector[1] + (1-x_per)*detector[0]
    return x

def circle_dist(detector,x_per,E):
    lam = E2lambda(E)
    theta = np.arcsin(lam/(2*d_SI111))
    inner_angle = np.pi - 2*theta
    
    x = pos_det(x_per,detector)
    D = np.linalg.norm(x)
    nu = np.arctan2(x[1],x[0])
    R = D*np.sin(np.pi/2-inner_angle)/np.sin(2*inner_angle)
    
    R_c = np.array([R*np.cos(nu-(np.pi/2-inner_angle)),R*np.sin(nu-(np.pi/2-inner_angle))])
    
    return R_c, R

# Error functions
# 1. Based on the distance to the circle center
def error_R(guess,x,y,detector,distribution,par):
    E = distribution(guess,par)
    R_c, R = circle_dist(detector,guess,E)
    dist_R = np.sqrt((x-R_c[0])**2 + (y-R_c[1])**2)
    return np.abs(dist_R-R)

# 2. Based on the energy and the angle, independent of the bragg circles
def error_theta(guess,x,y,detector,distribution,par):
    E = distribution(guess,par)
    v_1 = np.array([x,y])
    v_2 = pos_det(guess,detector) - v_1
    theta = np.arccos(np.dot(v_2,v_1)/(np.linalg.norm(v_1)*np.linalg.norm(v_2)))
    theta_bragg = (np.arcsin(E2lambda(E)/(2*d_SI111)))*2
    return np.abs(theta-theta_bragg)

# Find angles and energy in space
def points_to_angles(detector,pars,distribution=energy_det,mm=0.001):
    x_space, y_space = points_in_space(mm)
    angles = np.zeros((len(x_space),len(y_space)))
    energy_space = np.zeros((len(x_space),len(y_space)))    
    
    # Check if minimum of energy distribution is at the start, and the maximum at the end
    min_E = minimize(lambda x: distribution(x,pars),0,bounds=[(0,1)]).x[0]
    max_E = minimize(lambda x: -distribution(x,pars),0,bounds=[(0,1)]).x[0]
    max_E = np.round(max_E,5)
    if min_E != 0 or max_E != 1:
        print('Error, energy distribution not correct')
        return None, None, None, None
    
    for i in tqdm(range(len(x_space))):
        for j in range(len(y_space)):
            # Check if the point is in the area of interest
            if y_space[j] < -np.tan(2*np.pi/180)*x_space[i] or y_space[j] > np.tan(2*np.pi/180)*x_space[i]:
                angles[i,j] = np.nan
                energy_space[i,j] = np.nan
                continue
            
            # Find the correct energy
            per = minimize(lambda x: error_theta(x,x_space[i],y_space[j],detector,distribution,pars),0.5,bounds=[(0,1)]).x[0]
            
            # Check if minimization converged
            if error_theta(per,x_space[i],y_space[j],detector,distribution,pars) > 1e-5:
                angles[i,j] = np.nan  
                energy_space[i,j] = np.nan 
                continue
            
            # Calculate the angle
            E = distribution(per,pars)
            detector_ij = pos_det(per,detector)
            
            D_vec = np.array([detector_ij[0]-x_space[i],detector_ij[1]-y_space[j]])
            horizontal = np.array([1,0])
            theta_H = np.arccos(np.dot(D_vec,horizontal)/(np.linalg.norm(D_vec)*np.linalg.norm(horizontal)))
            theta_bragg = np.arcsin(E2lambda(E)/(2*d_SI111))
            angles[i,j] = theta_H-theta_bragg
            energy_space[i,j] = E
            
    # Limit the space to the area of interest
    angle_no_nan = np.copy(angles)
    angle_no_nan[np.isnan(angles)] = 0
    
    sum_y = np.sum(angle_no_nan,axis=0)
    sum_x = np.sum(angle_no_nan,axis=1)
    
    sum_y[sum_y==0] = np.nan
    sum_x[sum_x==0] = np.nan
    
    sum_y = pd.Series(sum_y)
    sum_x = pd.Series(sum_x)

    bottom = sum_y.first_valid_index()
    top = sum_y.last_valid_index()
    left = sum_x.first_valid_index()
    right = sum_x.last_valid_index()
    x_space = x_space[left:right]
    y_space = y_space[bottom:top]
    angles = angles[left:right]
    angles = angles[:,bottom:top]
    energy_space = energy_space[left:right]
    energy_space = energy_space[:,bottom:top]
    
    
    # Check for intersections
    for i in range(len(y_space)):
        energy = energy_space[:,i]
        first_idx = np.where(energy != np.nan)[0][0]
        last_idx = np.where(energy != np.nan)[0][-1]
        max_idx = np.where(energy == np.max(energy))[0]
        min_idx = np.where(energy == np.min(energy))[0]
        
        if first_idx != max_idx:
            print('Error, maximum energy not at the start')
            return None, None, None, None
        if last_idx != min_idx:
            print('Error, minimum energy not at the end')
            return None, None, None, None
    
    
    print('Angles in space calculated')    
    return x_space, y_space, angles, energy_space

#TODO
# Not fixed yet
def solutions(folder,x_space,y_space,angles,energies,spacing=0.01, offcut=0):
    angles = np.flip(angles.T,axis=0)
    angles = angles + offcut*np.pi/180
    dy = np.sin(angles)
    dx = np.cos(angles)
    X, Y = np.meshgrid(x_space,y_space)
    
    stream_points_x = np.arange(x_space[0],x_space[-1],spacing)

    stream_points_y = np.zeros_like(stream_points_x)
    stream_points = np.array([stream_points_x,stream_points_y]).T
    
    # Streamplot method
    plt.figure()
    
    plt.plot(x_space, np.tan(2*np.pi/180)*x_space,'k--',alpha=0.5)
    plt.plot(x_space, -np.tan(2*np.pi/180)*x_space,'k--',alpha=0.5)
    output = plt.streamplot(X,Y,dx,dy,density=10,start_points=stream_points,arrowsize=0)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(f'{folder}/streamplot.png')
    plt.close()
    output = output.lines.get_segments()
    output = split_lines(output)
    return output

def split_lines(output):
    k = 0
    line_collection = {}
    lines_x = []
    lines_y = []
    for i in range(len(output)):
        line = output[i]
        
        if i == 0:
            lines_x.append(line[0][0])
            lines_x.append(line[1][0])
            lines_y.append(line[0][1])
            lines_y.append(line[1][1])
            continue
        if i == len(output)-1:
            lines_x = np.array(lines_x)
            lines_y = np.array(lines_y)
            line_collection[k] = np.array([lines_x,lines_y])
            continue
        
        if line[0][0] == lines_x[-1] and line[0][1] == lines_y[-1]:
            lines_x.append(line[1][0])
            lines_y.append(line[1][1])
        else:
            lines_x = np.array(lines_x)
            lines_y = np.array(lines_y)
            line_collection[k] = np.array([lines_x,lines_y])
            k += 1
            lines_x = []
            lines_y = []
            lines_x.append(line[0][0])
            lines_x.append(line[1][0])
            lines_y.append(line[0][1])
            lines_y.append(line[1][1])
    
    return line_collection

def circle_fit(line):
    line = line.T
    R_x,R_y, R, sigma = cf.least_squares_circle(line)
    return R_x, R_y, R, sigma

def mcstas_analyzers(folder,mm, offcut,output, mos=6):
    def distance(x1,y1,x2,y2):
        return np.sqrt((x1-x2)**2 + (y1-y2)**2)
    
    string = ""
    height_deg = 2.5
    for key in output.keys():
        line = output[key]
        R_x, R_y, R, sigma = circle_fit(line)
        
        if R > 5:
            R = 0
        
        # Define points within +- 2 degrees off of line [0,0] to [5,0]
        
        x_c, y_c = circle_in_area(R_x,R_y,R,angle=2)
        dist = x_0(R_x,R_y,R)
        angle = np.pi/2 - np.arctan2(R_y,dist-R_x)
        
        angle %= np.pi
        
        width = distance(x_c[0],y_c[0],x_c[-1],y_c[-1]) 
        height = 2*np.tan(height_deg*np.pi/360)*line[0,0]

        d = {"Number": key, "width": width, "height": height, "mm": mm/2, "radius": R, "reflection": ref, "offcut": offcut, "mosaicity": mos, "Dist": dist, "Angle": angle*180/np.pi}
        analyzer_text = ""
        with open(f"{main_folder}/WARP_Single_Analyzer.txt",'r') as f:
            text = f.readlines()
            for line in text:
                if "{" in line:
                    s = line.find("{")
                    f = line.find("}")
                    text_key = line[s+1:f]
                    analyzer_text += line[:s] + str(d[text_key]) + line[f+1:]
                    
                else:
                    analyzer_text += line  
        string += analyzer_text
        string += "EXTEND %{\n"
        string += "    if(SCATTERED){scat=1;}\n"
        string += "%}\n"
        string += "\n\n"
    return string

def mcstas_instr(folder,d_pars,E_pars, analyzer_text):
    E_min, E_max, *args = E_pars
    L0, L1, Ld, angle = d_pars
    
    E_0 = 0.5*(E_max+E_min)
    d_E = 0.5*(E_max-E_min)
    
    # Split analyzer_text
    analyzers = analyzer_text.split("\n")
    
    with open(f"{main_folder}/WARP_Instr_Template.instr",'r') as f:
        text = f.readlines()
        text[28] = text[28][:-2] + f"{E_0},\n"
        text[29] = text[29][:-2] + f"{d_E},\n"
        text[33] = text[33][:-2] + f"{L0},\n"
        text[34] = text[34][:-2] + f"{L1},\n"
        text[35] = text[35][:-2] + f"{Ld},\n"
        text[36] = text[36][:-2] + f"{angle},\n"
        string = '"x bins={} energy limits=[{} {}] bins={}"'.format(Ld//0.001,E_min*0.9,E_max*1.1,(E_max-E_min)//0.001)
        text[112] = text[112][:-2] + f"{string},\n"
        
        # Add analyzers at line 121
        for i in range(len(analyzers)):
            text.insert(90+i,analyzers[i]+ "\n")
        
    
    with open(f"{folder}/WARP.instr",'w') as f:
        f.writelines(text)
    
    shutil.copy(f"{main_folder}/Monochromator_bent.comp",f"{folder}/Monochromator_bent.comp")
    
    print('McStas Instrument file created')
    return True


# Plotting functions
def plot_energy_dist(folder,E_pars):
    E_max, E_min, dE_0, dE_1, ddE_0 = E_pars
    x = np.linspace(0,1,1000)
    E = energy_det(x,E_pars)
    plt.figure()
    plt.plot(x,E)
    plt.xlabel('x')
    plt.ylabel('Energy')
    plt.title('Energy distribution')
    plt.savefig(f'{folder}/energy_dist.png')
    plt.close()

def points_2D_map(folder,x,y,angle, name='angles'):    
    angle = np.flip(angle.T,axis=0)
    if name == "angles":
        angle = angle*180/np.pi
    plt.figure()
    plt.imshow(angle,extent=(x[0],x[-1],y[0],y[-1]),aspect='auto')

    plt.plot(x, np.tan(2*np.pi/180)*x,'k--',alpha=0.5)
    plt.plot(x, -np.tan(2*np.pi/180)*x,'k--',alpha=0.5)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Angles in space')
    plt.savefig(f'{folder}/{name}.png')
    plt.close()

def reflection_map(folder,x,y,reflection,energy, d_pars,E_pars, detector):
    E_max, E_min, dE_0, dE_1, ddE_0 = E_pars
    pos = np.linspace(0,1,10,endpoint=True)
    d_linspace = []
    for p in pos:
        d_linspace.append(pos_det(p,detector))
    E_linspace = energy_det(pos,E_pars)
    L0, L1, Ld, angle = d_pars
    xs = []
    ys = []
    es = []
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(E_linspace)):
                if math.isclose(energy[i,j],E_linspace[k],abs_tol=1e-2):
                    xs.append(x[i])
                    ys.append(y[j])
                    es.append(E_linspace[k])
    arg_sort = np.argsort(es)
    xs = np.array(xs)[arg_sort]
    ys = np.array(ys)[arg_sort]
    es = np.array(es)[arg_sort]                
    
    # Energy color map
    # High energy is red, low energy is blue
    colors = (E_linspace - E_min)/(E_max-E_min)
    colors = colors[::-1]
    colors = plt.cm.coolwarm(colors)
    
    # Plot the reflection map, including sample, and detector
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(0,0,'ok',label='Sample')
    ax.plot([detector[0][0],detector[1][0]],[detector[0][1],detector[1][1]],'k-',label='Detector',linewidth=2.5)
    idx_used = []
    for i in tqdm(range(len(xs))):
        # Plot energy
        idx = np.where(E_linspace == es[i])[0][0]
        if idx not in idx_used:
            ax.plot([0,xs[i]],[0,ys[i]],color=colors[idx],linewidth=0.5,label=f"Energy: {es[i]:.2f}")
            ax.plot([xs[i],d_linspace[idx][0]],[ys[i],d_linspace[idx][1]],color=colors[idx],linewidth=1)
            idx_used.append(idx)
        else:
            ax.plot([0,xs[i]],[0,ys[i]],color=colors[idx],linewidth=0.5)
            ax.plot([xs[i],d_linspace[idx][0]],[ys[i],d_linspace[idx][1]],color=colors[idx],linewidth=1)
    
    # Make scale equal for x and y
    ax.set_aspect('equal')
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Reflection map')
    ax.legend()
    plt.savefig(f'{folder}/reflection_map.png')
    plt.close()

def curvature_map(folder,output):
    R_c = []
    S_c = []
    for key in output.keys():
        R_x, R_y, R, sigma = circle_fit(output[key])
        R_c.append(R)
        S_c.append(sigma)
    fig, ax = plt.subplots()
    
    ax.plot(R_c, color='blue')
    ax.set_ylabel('Radius of Curvature / m', color='blue')
    ax.set_xlabel('Circle number')
    ax.set_yscale('log')
    ax.hlines(5,0,len(R_c),color='black',linestyles='dashed')
    
    ax2 = ax.twinx()
    ax2.set_ylabel('Error on Curvature', color='red')
    ax2.plot(S_c, color='red')
    ax2.set_yscale('log')
    plt.savefig(f'{folder}/Radius.png')
    plt.close()

def plot_circle(line):
    R_x, R_y, R, sigma = circle_fit(line)
    theta = np.linspace(0,2*np.pi,100)
    x = R_x + R*np.cos(theta)
    y = R_y + R*np.sin(theta)
    plt.plot(x,y)
    plt.plot(R_x,R_y,'o')
    plt.plot(line[0],line[1])
    plt.show()

def plot_R_c(folder,output):
    R_x = []
    R_y = []
    for key in output.keys():
        R_x_temp, R_y_temp, R, sigma = circle_fit(output[key])
        R_x.append(R_x_temp)
        R_y.append(R_y_temp)
    plt.figure()
    plt.plot(R_x,R_y)
    plt.plot(0,0,'x')
    x = np.linspace(0,5,100)
    plt.plot(x,np.tan(2*np.pi/180)*x,'k--',alpha=0.5)
    plt.plot(x,-np.tan(2*np.pi/180)*x,'k--',alpha=0.5)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Circle centers')
    plt.savefig(f'{folder}/R_c.png')
    plt.close()

# Main function
def start(d_pars,E_pars, mm):

    L0, L1, Ld, angle = d_pars

    # Distribute the detector
    det0 = np.array([L0,L1])
    d_det = np.array([np.cos(angle),np.sin(angle)])*Ld/2
    detector = np.array([det0-d_det,det0+d_det])
    
    return detector

def save_data(folder, x, y, angles, energies):
    
    if x is None:
        return False
    
    # Make folder

    os.makedirs(folder,exist_ok=True)
    
    # Save plots
    points_2D_map(folder, x,y,angles, name='angles')
    points_2D_map(folder, x,y,energies, name=f'energies')
    
    # Save data
    np.save(f'{folder}/x.npy',x)
    np.save(f'{folder}/y.npy',y)
    np.save(f'{folder}/angles.npy',angles)
    np.save(f'{folder}/energies.npy',energies)
    
    return True

def import_data(folder,n):
    x = np.load(f'{folder}/x.npy')
    y = np.load(f'{folder}/y.npy')
    angles = np.load(f'{folder}/angles.npy')
    energies = np.load(f'{folder}/energies.npy')
    
    with open(f"{main_folder}/parameter.txt", "r") as f:
        lines = f.readlines()
        line = lines[n+1]
        line = line.split(";")
        E_pars = [float(i) for i in line[6:11]]
        d_pars = [float(i) for i in line[2:6]]
        mm = float(line[11])
    
    return x, y, angles, energies, d_pars, E_pars, mm

def parameter_file(n, d_pars,E_pars,mm, convergence):
    L0, L1, Ld, angle = d_pars
    E_max, E_min, b, rot, offset = E_pars
    
    if convergence:
        conv = "Yes"
    else:
        conv = "No"
    
    with open(f"{main_folder}/parameter.txt", "a") as f:
        date_time = datetime.datetime.now().strftime("%Y-%m-%d")
        f.write(f"{n}; {date_time}; {L0}; {L1}; {Ld}; {angle}; {E_max}; {E_min}; {b}; {rot}; {offset}; {mm}; {conv}\n")
    
    return n

def run(folder,d_pars,E_pars,mm, offcut, spacing, Load=0, mcstas=True):
    # Make detector distribution
    detector = start(d_pars,E_pars,mm)
    # Make folder
    
    # Start the process
    if not Load:
        n, folder = gen_folder(folder)
        # Generate the grid, angles and energies
        x,y,angles,energies = points_to_angles(detector,E_pars,mm=mm)
        if x is not None:
            save_data(folder, x, y, angles, energies)
            n = parameter_file(n, d_pars,E_pars,mm, True)
            
        else:
            n = parameter_file(-1, d_pars,E_pars,mm, False)
            return False
    else:
        folder = f'{main_folder}/Data/{Load}'
        n = Load
        x,y,angles,energies, d_pars, E_pars, mm = import_data(folder, n)
        print("Done Loading")
    print(f"Folder nr. {n}")
    # Generate Solutions
    plot_energy_dist(folder,E_pars)
    reflection_map(folder,x,y,ref,energies, d_pars,E_pars, detector)
    output = solutions(folder,x,y,angles,energies,spacing,offcut=offcut)
    curvature_map(folder,output)
    plot_R_c(folder,output)
    print("Plotted all")
    # Generate McStas files
    analyzer_text = mcstas_analyzers(folder,mm, offcut,output, mos=20)
    mcstas_instr(folder,d_pars,E_pars,analyzer_text)
    folder = f'{main_folder}/Data/{n}'
    string = f"C:/mcstas-3.4/bin/mcrun -c -n1e6 {folder}/WARP.instr -d {folder}/Test focus=4"
    
    
    if mcstas:
        if os.path.exists(f"{folder}/Test/mccode.sim"):
            print("Simulation already run")
        else:
            print("Running: McStas simulation")
            call(string,shell=True)
        
        # Plot the results
        mcstas_x_E_plot(folder)
        mcstas_resolution_function(folder,d_pars,E_pars)
        
        
    return n


# Constants
main_folder = "Marmot_Opt"
ref = '"Si111"'     # Reflection
d_SI111 = 3.13536   # d-spacing of Si111
mos = 6             # Mosaicity
spacing = 0.01  # Spacing between analyzers


# Marmot Distances
top = np.array([89.94,222.9])
bottom = np.array([144.2,216.87])
L0 = ((top[0]+bottom[0])/2)/100 # In meters
L1 = ((top[1]+bottom[1])/2)/100 # In meters
Ld = (np.linalg.norm(top-bottom))/100 # In meters
d_angle = -(-np.arctan2(top[1]-bottom[1],top[0]-bottom[0])*180/np.pi % 180)
d_pars = (L0,L1,Ld,d_angle)
print(d_pars)

# Marmot Energy
E_min = 3.95
E_max = 6.71
dE_0 = 1
dE_1 = 4
ddE_0 = 3
E_pars = (E_min,E_max,dE_0,dE_1,ddE_0)

# WARP Distance parameters
#L0 = 2
#L1 = 3.5
#Ld = 1
#angle = -5
#d_pars = (L0,L1,Ld,angle)


# WARP Energy parameters
#E_0 = 2.5
#E_1 = 5
#dE_0 = 0.2 # Has to be positive
#dE_1 = 4 # Has to be positive
#ddE_0 = 0.1
#E_pars = (E_0,E_1,dE_0,dE_1,ddE_0)

# Analyzer parameters
mm = 0.001      # Spacing between points in generated grid
offcut = 0      # Offcut in degrees
main_folder = "Marmot"

n = run(main_folder,d_pars,E_pars,mm, offcut, spacing, Load=0, mcstas=False)


"""
if spacing < mm:
    print("Spacing less than mm, need finer spacing in grid")
    raise ValueError

dE_0_ar = np.linspace(0.1,2,20)

os.makedirs(main_folder,exist_ok=True)
shutil.copy(f"Marmot/WARP_Single_Analyzer.txt",f"{main_folder}/WARP_Single_Analyzer.txt")
shutil.copy(f"Marmot/WARP_Instr_Template.instr",f"{main_folder}/WARP_Instr_Template.instr")
shutil.copy(f"Marmot/Monochromator_bent.comp",f"{main_folder}/Monochromator_bent.comp")
with open(f"{main_folder}/parameter.txt", "w") as f:
    f.write("# Which parameters have been run through\n# Run; Time; d_pars; E_pars; mm; Converged?\n")
for dE_0 in dE_0_ar:
    E_pars = (E_min,E_max,dE_0,dE_1,ddE_0)
    n = run(main_folder,d_pars,E_pars,mm, offcut, spacing, Load=0, mcstas=True)
"""


# TODO
# 7. Make a array that records time for each slice in space
# 8. Make parameter space (L1, L0, Ld, E_max, E_min, b, rot, offset)
# 9. Scan over parameter space
# 10. Different offcuts for different lines

