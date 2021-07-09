from math import floor,sqrt,exp
import numpy as np
import sys
sys.path.append('/Users/chrispdavies/Documents/University/University of \
Bath/Hill_Profiling/GivenPrograms/Chris_Davies_hill-profiling_Project_2017')
from CEDA_dtm_5m_01 import DTM

dtm = DTM()


def interp_surf(x,y,H) :
    # 2x2 matrix of heights at four grid points, H
    # x and y are length 2 vectors describing the grid points
    # interp_surf returns a vector of the coefficents of the 2D function:
    # axy + bx + cy + d = z
    # which is the surface gained by bilinear interpolation
    f = ( x[1] - x[0] )*( y[1] - y[0] )
    a = ( H[0,0] + H[1,1] - H[0,1] - H[1,0] )/f
    b = ( H[0,1]*y[0] + H[1,0]*y[1] - H[0,0]*y[1] - H[1,1]*y[0] )/f
    c = ( x[1]*( H[0,1] - H[0,0] ) + x[0]*( H[1,0] - H[1,1] ) )/f
    d = ( x[1]*( H[0,0]*y[1] - H[0,1]*y[0] ) + x[0]*( H[1,1]*y[0] - \
    H[1,0]*y[1] ) )/f
    return (a,b,c,d)

def terrain(x,y,f) :
    # given a point in the x-y plane and surface described by f, 
    # terrain returns the height of the surface at that point
    return f[0]*x*y + f[1]*x + f[2]*y + f[3]

def line_surf_inter(p,q,f,x,y) :
    # p and q are points in 3D space
    # f is a vector of the coefficents of the function axy + bx + cy + d = z
    # which describes a surface. 
    # x and y describe the coordinates of the four grid points
    # line_surf_inter determines whether the line between the points 
    # intersects the surface within the four grid points

    # is one of the points on the plane and within the grid points?
    fp = f[0]*p[0]*p[1] + f[1]*p[0] +f[2]*p[1] + f[3]
    fq = f[0]*q[0]*q[1] + f[1]*q[0] +f[2]*q[1] + f[3]

    if  abs(fp-p[2]) < 1e-6 :
        if ((x[0] <= p[0] <= x[1]) and (y[0] <= p[1] <= y[1])) :
            return True
    
    if abs(fq-q[2]) < 1e-6 :
        if ((x[0] <= q[0] <= x[1]) and (y[0] <= q[1] <= y[1])) :
            return True

    a = q-p 
    u = f[0]*a[0]*a[1]
    v = f[0]*( p[0]*a[1] + p[1]*a[0] ) + f[1]*a[0] + f[2]*a[1] - a[2]
    w = f[0]*p[0]*p[1] + f[1]*p[0] + f[2]*p[1] + f[3] - p[2]
    # line between p and q is s = p + (q-p)*t  for paramter t
    # substituting the line equation for the x,y, and z values in the 
    # surface equation yields a quadratic in t: ut^2 + vt + w = 0
    if v**2 < 4*u*w :
        return False # no real t solution => line does not intersect surface
    elif abs(u) < 1e-6 :
        if abs(v**2)  < 1e-6 :
            return False
        else:
            t = -w/v
            if 0 <= t <= 1 :
                s = p + a*t
                return ((x[0] <= s[0] <= x[1]) and (y[0] <= s[1] <= y[1]))
    else :
        t1 = (-v + (v**2 - 4*u*w)**.5)/(2*u)
        t2 = (-v + (v**2 - 4*u*w)**.5)/(2*u)
        tol = 5e-4 
        # is one of the intersection points between p and q?
        if 0 <= t1 <= 1 :
            s = p + a*t1 # intersection point
            # is the intersection point between the grid points?
            if (x[0]-tol <= s[0] <= x[1]+tol) and (y[0]-tol <= s[1] <= y[1]+tol) :
                return True 
        elif 0 <= t2 <= 1 :
            r = p + a*t2 # intersection point

            if (x[0]-tol <= r[0] <= x[1]+tol) and (y[0]-tol <= r[1] <= y[1]+tol) :
                return True
        
        else :
            return False

def grid_squ(p,q,x,y) :
    # grid described by vectors x and y 
    # points in 3D space p and q
    # grid_sq returns a nx2 array describing the lower left grid points  of the 
    # grid squares through which the line between p and q passes 
    
    xd = x[1]-x[0]
    yd = y[1]-y[0]
    xl = x[0]
    yl = y[0]
    a = q - p
    
    # Establish minimum and maximum x coordinates and grid lines
    pq_x = np.array([p[0],q[0]])
    pqmin_x = np.amin(pq_x) 
    pqmax_x = np.amax(pq_x)
    gmin_x = x[floor((pqmin_x-xl)/xd)]
    if abs(gmin_x+xd-pqmin_x) < 1e-8 :
        gmin_x += xd
    gmax_x = x[floor((pqmax_x-xl)/xd)]
    if abs(gmax_x+xd-pqmax_x) < 1e-8 :
        gmax_x += xd
    c_x = gmin_x + xd - pqmin_x 
    d_x = pqmax_x - gmax_x 

    # Establish minimum and maximum y coordinates and grid lines
    pq_y = np.array([p[1],q[1]])
    pqmin_y = np.amin(pq_y)
    pqmax_y = np.amax(pq_y)
    gmin_y = y[floor((pqmin_y-yl)/yd)]
    if abs(gmin_y+yd-pqmin_y) < 1e-8 :
        gmin_y += yd
    gmax_y = y[floor((pqmax_y-yl)/yd)]
    if abs(gmax_y+yd-pqmax_y) < 1e-8 :
        gmax_y += yd
    c_y = gmin_y + yd - pqmin_y
    d_y = pqmax_y - gmax_y 

    n_x = int(round(1 + (abs(a[0]) - c_x - d_x)/xd)) #number of intersections 
    #with x grid lines
    n_y = int(round(1 + (abs(a[1]) - c_y -d_y)/yd)) #number of intersections
    #with y grid lines

    n = int(n_x + n_y + 1) # maximum number of grid squares

    c = np.zeros((n,2))


    if abs(a[0]) < 1e-8 :
        c[:,0] = gmin_x 
        c[:,1] = np.arange(gmin_y,gmin_y+(n_y+.5)*yd,yd)  
    elif abs(a[1]) < 1e-8 :
        c[:,1] = gmin_y 
        c[:,0] = np.arange(gmin_x,gmin_x+(n_x+.5)*xd,xd)
    else :

        if abs(d_x) < 1e-8 :
            gmax_x += -xd
            n_x += -1
    
        if abs(d_y) < 1e-8 :
            gmax_y += -yd
            n_y += -1

        n = int(n_x + n_y + 1) # maximum number of grid squares

        c = np.zeros((n,2))

        x0 = pqmin_x 
        t = (x0 - p[0])/a[0]
        y0 = p[1]+a[1]*t
        grid_y0 = y[floor((y0-yl)/yd)]
        if abs(grid_y0+yd-y0) < 1e-8 :
            grid_y0 += yd 

        t = 1.0 - t 
        y2 = p[1]+a[1]*t 
        grid_y2 = y[floor((y2-yl)/yd)]
        if abs(grid_y2+yd-y2) < 1e-8 :
            grid_y2 += yd

        j = 0

        if abs(d_y) < 1e-8 : 
            if grid_y2 > gmin_y :
                grid_y2 += -yd
            else:
                grid_y0 += -yd 

        for i in range(0,n,1) :
            c[i,0] = gmin_x 
            c[i,1] = grid_y0

            if gmin_x == gmax_x and grid_y0 == grid_y2 :
                break

            t = 1.0*(gmin_x + xd - p[0])/a[0]
            yg = p[1]+a[1]*t


            if abs(yg - grid_y0) < 1e-8 :
                grid_y0 += -yd
                gmin_x += xd
                j += 1
            elif abs(yg - (grid_y0 + yd)) < 1e-8 :
                grid_y0 += yd
                gmin_x += xd
                j += 1
            elif yg < grid_y0 :
                grid_y0 += -yd
            elif grid_y0 < yg < grid_y0 + yd :
                gmin_x += xd
            elif yg > grid_y0 + yd:
                grid_y0 += yd 
            
        k =  np.arange(n-j,n-0.5,1)
        c = np.delete(c,k,0)
    return c 

def intervis(p,q,x,y,H) :
    # grid described by vectors x and y
    # H matrix of heights of grid points
    # points p and q
    # intervis returns True if the points p and q are intervisible, 
    # False otherwise 
    c = grid_squ(p,q,x,y) 
    if c.shape[0] < 1 :
        print(q)
    n = c.shape[0]
    xd = x[1]-x[0]
    yd = y[1]-y[0]

    pq_z = np.array([p[2],q[2]])
    pqmin_z = np.amin(pq_z)
    j = int((c[0,0]-x[0])/xd)
    k = int((c[0,1]-y[0])/yd)
    l = int((c[n-1,0]-x[0])/xd)
    m = int((c[n-1,1]-y[0])/yd) 
    km = np.array([k,m])
    k = np.amin(km)
    m = np.amax(km)
    H1 = H[j:l+2,k:m+2]
    H_max = np.amax(H1)

    if pqmin_z > H_max :
        return True 
    else :

        for i in range(0,n,1) :
            j = int((c[i,0]-x[0])/xd)
            k = int((c[i,1]-y[0])/yd)
            x1 = np.array([c[i,0],c[i,0]+xd])
            y1 = np.array([c[i,1],c[i,1]+yd])
            H1  = H[j:j+2,k:k+2]
            if H1.shape != (2,2) :
                print(i)
            f = interp_surf(x1,y1,H1)
            if line_surf_inter(p,q,f,x1,y1) :
                #print(c[i,0],c[i,1],i)
                return False
        return True 

def hor_angle(p,a,x,y,H) :
    # 3D point p, direction in the x-y plane a,
    # grid described by vectors, x and y, matrix of heights at grid points H
    nx = x.shape[0]
    ny = y.shape[0]
    axy = (a[0]**2 + a[1]**2)**.5
    phi = np.arctan2(a[1],a[0])
    if abs(a[0]) < 1e-8 :
        if a[1] < 0 :
            d = p[1]-y[1]
        elif a[1] > 0 :
            d = y[ny-2]-p[1]
    elif abs(a[1]) < 1e-8  :
        if a[0] > 0 :
            d = x[nx-2]-p[0]
        else :
            d = p[0]-x[1]
    elif -np.pi < phi < -np.pi/2 :
        d = axy*np.amin([abs((p[0]-x[1])/a[0]),abs((p[1]-y[1])/a[1])])
    elif -np.pi/2 < phi < 0 :
        d = axy*np.amin([(x[nx-2]-p[0])/a[0],abs((p[1]-y[1])/a[1])])
    elif 0 < phi < np.pi/2 :
        d = axy*np.amin([(x[nx-2]-p[0])/a[0],(y[ny-2]-p[1])/a[1]])
    elif np.pi/2 < phi < np.pi :
        d = axy*np.amin([abs((p[0]-x[1])/a[0]),(y[ny-2]-p[1])/a[1]])
      
    theta = np.pi/2
    alpha = theta/4
    theta2 = 20
    a1 = a/axy 
    if abs(a1[0]) < 1e-8 :
        a1[0]=0
    if abs(a1[1]) < 1e-8 :
        a1[1]=0
    q = np.zeros(3)
    while theta >= -np.pi/2 :  
        q[2] = p[2] + d*np.sin(theta)
        q[0] = p[0] + d*np.cos(theta)*a1[0]
        q[1] = p[1] + d*np.cos(theta)*a1[1]
        if intervis(p,q,x,y,H) :
            theta2 = theta
            theta += - alpha 
        else :
            break 
    if theta2 == 20 :
        return theta
    theta1 = theta2 - alpha/2
    while d*np.sin(alpha) > 0.5 :
        q[2] = p[2] + d*np.sin(theta1)
        q[0] = p[0] + d*np.cos(theta1)*a1[0]
        q[1] = p[1] + d*np.cos(theta1)*a1[1]
        if intervis(p,q,x,y,H) :
            theta2 = theta1
            alpha *= 0.5
            theta1 = theta2 - alpha/2
        else :
            theta = theta1
            alpha *= 0.5
            theta1 = theta2 - alpha/2
    return theta1 


def horizon(p,a,theta,n,x,y,H) :
    # point p, direction in x-y plane a, angle theta in radians between 0 and pi
    # grid described by vectors x and y, matrix of heights at grid points
    
    b = np.arctan2(a[1],a[0])
    b1 = b - theta/2
    b2 = b + theta/2
    phi = np.arange(b1,b2+theta/(2*n),theta/(n-1))
    ax = np.cos(phi)
    ay = np.sin(phi)
    alpha = np.zeros(phi.shape)

    for i in range(0,phi.shape[0],1) :
        a1 = np.array([ax[i],ay[i]])
        alpha[i] = hor_angle(p,a1,x,y,H)
    return alpha,phi

def intervis2(p,q) :
    # given two points and a surface described by the function elevation, 
    # intervis2 determines whether the two points are intervisible
    a = q-p 
    m = p + a*0.5
    n = np.round(m,0)
    n = n.astype(int)
    if np.all(abs(n-p)<1e-8) or np.all(abs(n-q)<1e-8) :
        return True
    elif sqrt(a[0]**2 + a[1]**2) <= 1.0 :
        return True
    #elif n[2] <= dtm.height(n[0],n[1]) :
    #    return False
    elif n[2] <= elevation(n[0],n[1]) :
        return False 
    elif intervis2(p,m) and intervis2(m,q) :
        return True 
    else :
        return False

def elevation(x,y):
  # four hills of heights h0,h1,h2,h3 at (x0,y0),(x1,y1),(x2,y2) and (x3,y3)
    x0,y0=92.0,23.0
    x1,y1=86.0,76.0
    x2,y2=25.0,83.0
    x3,y3=150.0,122.0
    h0,h1,h2,h3=60.0,50.0,65.0,80.3
    w0,w1,w2,w3=30.0,30.0,30.0,40.0 # width parameters 
    return h0*np.exp(-(np.hypot(np.subtract(x,x0),np.subtract(y,y0))/w0)**2) + \
        h1*np.exp(-(np.hypot(np.subtract(x,x1),np.subtract(y,y1))/w1)**2) + \
        h2*np.exp(-(np.hypot(np.subtract(x,x2),np.subtract(y,y2))/w2)**2) + \
        h3*np.exp(-(np.hypot(np.subtract(x,x3),np.subtract(y,y3))/w3)**2)

def hor_angle2(p,a,d) :
    # Given a point in 3D space, p, a direction, a, and a distance, d, 
    # hor_angle2 gives the angle to the horizon in the direction a
    # terrain described by elevation

    axy = sqrt(a[0]**2 + a[1]**2)
    a1 = a/axy
    theta = np.pi/2
    alpha = theta/4
    theta2 = 20
    a1 = a/axy 
    q = np.zeros(3)
    while theta >= -np.pi/2 :  
        q[2] = p[2] + d*np.sin(theta)
        q[0] = p[0] + d*np.cos(theta)*a1[0]
        q[1] = p[1] + d*np.cos(theta)*a1[1]
        if intervis2(p,q) :
            theta2 = theta
            theta += - alpha 
        else :
            break 
    if theta2 == 20 :
        return theta
    theta1 = theta2 - alpha/2
    while d*np.sin(alpha) > 0.5 :
        q[2] = p[2] + d*np.sin(theta1)
        q[0] = p[0] + d*np.cos(theta1)*a1[0]
        q[1] = p[1] + d*np.cos(theta1)*a1[1]
        if intervis2(p,q) :
            theta2 = theta1
            alpha *= 0.5
            theta1 = theta2 - alpha/2
        else :
            theta = theta1
            alpha *= 0.5
            theta1 = theta2 - alpha/2
    return theta1 

def horizon2(p,a,theta,n,d) :
    # point p, direction in x-y plane a, angle theta in radians between 0 and pi
    # terrain described by elevation, horizon2 returns a vector of length n+1 
    # of angles describing the horizon theta/2 either side of the direction a
    
    b = np.arctan2(a[1],a[0])
    b1 = b - theta/2
    b2 = b + theta/2
    phi = np.arange(b1,b2+theta/(2*n),theta/(n-1))
    ax = np.cos(phi)
    ay = np.sin(phi)
    alpha = np.zeros(phi.shape)

    for i in range(0,phi.shape[0],1) :
        a1 = np.array([ax[i],ay[i]])
        alpha[i] = hor_angle2(p,a1,d)
    return alpha,phi

def hill_shape(vertical,horizontal) :
    # hill_shape returns a vector of characteristics to describe the shape 
    # of the plot generated by horizon
    n = vertical.size 
    if n != horizontal.size :
        print('Vector sizes do not match')

    v_min = np.amin(vertical)
    v_max = np.amax(vertical)
    v_av = np.average(vertical)
    apparent_height = v_max - v_av 

    v = vertical - v_av
    v *= 1/apparent_height
    j = np.argmax(vertical) 
    v_max_h = horizontal[j]

    w = np.zeros(7)

    k = 0
    i0 = 0
    i1 = 0
    for i in range(0,j+1) :
        if v[j-i] <= 0.75 :
            if k == 0 :
                w[2] = (np.tan(v_max_h-horizontal[j-i]))/np.tan(apparent_height)
                i0 = j-i 
                k += 1
        if v[j-i] <= 0.5 :
            if k == 1 :
                w[1] = (np.tan(v_max_h-horizontal[j-i]))/np.tan(apparent_height)
                i0 = j-i
                k += 1
                break
        #if v[j-i] <= 0.25 :
        #    if k == 2 :
        #        w[1] = (np.tan(v_max_h-horizontal[j-i]))/np.tan(apparent_height)
        #        i0 = j-i
        #        k += 1
        #        break 
    
    if w[1] == 0 :
        w[1] = np.tan(horizontal[j]-horizontal[0])/np.tan(apparent_height)
    if w[2] == 0 :
        w[2] = np.tan(horizontal[j]-horizontal[0])/np.tan(apparent_height) 


    k = 0
    for i in range(j,n-1) :
        if v[i] <= 0.75 :
            if k == 0 :
                w[3] = (np.tan(horizontal[i]-v_max_h))/np.tan(apparent_height)
                i1 = i
                k += 1
        if v[i] <= 0.5 :
            if k == 1 :
                w[4] = (np.tan(horizontal[i]-v_max_h))/np.tan(apparent_height)
                i1 = i
                k += 1
                break
        #if v[i] <= 0.25 :
        #    if k == 2 :
        #        w[6] = (np.tan(horizontal[i]-v_max_h))/np.tan(apparent_height)
        #        i1 = i
        #        k += 1
        #        break

    if w[3] == 0 :
        w[3] = np.tan(horizontal[n-1]-horizontal[j])/np.tan(apparent_height)
    if w[4] == 0 :
        w[4] = np.tan(horizontal[n-1]-horizontal[j])/np.tan(apparent_height)

    gradients = np.zeros(n-1)
    for i in range(0,n-1) :
        gradients[i] = (np.tan(vertical[i+1]-vertical[i]))/(np.tan(horizontal[i+1]\
        -horizontal[i]))
    
    tol = 0.01
    k = 0
    r = np.zeros(n)
    s = np.zeros(n)
    if gradients[0] > tol :
        g0 = 1  
    elif gradients[0] < -tol :
        g0 = -1   
    else :
        g0 = 0
    r[0] = g0
    j = 0
    for i in range(1,n-1) :
        if gradients[i] > tol :
            g1 = 1  
        elif gradients[i] < -tol :
            g1 = -1   
        else :
            g1 = 0
        if g1 != g0 :
            d_h = (0.5*(horizontal[i]+horizontal[i+1]) - \
            0.5*(horizontal[j]+horizontal[j+1]))/apparent_height
            d_v = 0.5*(v[i]+v[i+1]) - 0.5*(v[j]+v[j+1])
            s[k] = sqrt(d_h**2 + d_v**2)
            j = i
            k += 1
            r[k] = g1
        g0 = g1
    d_h = (0.5*(horizontal[n-2]+horizontal[n-1]) - \
    0.5*(horizontal[j]+horizontal[j+1]))/apparent_height
    d_v = 0.5*(v[n-2]+v[n-1]) - 0.5*(v[j]+v[j+1])
    s[k] = sqrt(d_h**2 + d_v**2)

    for i in range(0,k+1) :
        if r[i] == 1 and s[i] > 0.3 :
            if r[i+1] == 0 :
                if r[i+2] == -1 and s[i+2] > 0.3 :
                    w[0] += 1

    j = np.argmax(vertical) 
    for i in range(1,j) :
        if gradients[j-i] > 0.05 :
            #w[7] = (-np.tan(horizontal[j-i]) - 0.5*(np.tan(horizontal[j-i+1])-\
            #np.tan(horizontal[j-i])) + np.tan(v_max_h))/np.tan(apparent_height) 
            w[5] = np.tan(-horizontal[j-i] - 0.5*(horizontal[j-i+1]-\
            horizontal[j-i]) + v_max_h)/np.tan(apparent_height) 
            break 
    for i in range(j,n-2) :
        if gradients[i] < -0.05 :
            #w[8] = (np.tan(horizontal[i]) + 0.5*(np.tan(horizontal[i+1])-\
            #np.tan(horizontal[i])) - np.tan(v_max_h))/np.tan(apparent_height)
            w[6] = np.tan(horizontal[i] + 0.5*(horizontal[i+1]-\
            horizontal[i]) - v_max_h)/np.tan(apparent_height)
            break 

    if w[5] == 0 :
        w[5] = np.tan(horizontal[j]-horizontal[0])/np.tan(apparent_height)
    if w[6] == 0 :
        w[6] = np.tan(horizontal[n-1]-horizontal[j])/np.tan(apparent_height)
    
    return w 

def piecewise_linear(x,c,e) :
    # given vectors c and e, piecewise_linear returns the value of the 
    # piecewise linear function described by c and e at the value x
    n = e.size 
    if n != c.size :
        print('Vector sizes do not match')
    x_l = np.amin(e)
    i = np.argmin(e)
    y_l = c[i]
    x_u = np.amax(e)
    i = np.argmax(e)
    y_u = c[i]
    for i in range(0,n) :
        if e[i] < x :
            if e[i] > x_l :
                if abs(e[i]-x) > 1e-6 :
                    x_l = e[i]
                    y_l = c[i]
        if e[i] > x :
            if e[i] < x_u :
                if abs(e[i]-x) > 1e-6 :
                    x_u = e[i]
                    y_u = c[i]
    if abs(x_u-x_l) < 1e-6 :
        y = y_l
    else :
        t = (x-x_l)/(x_u-x_l)
        y = y_l + t*(y_u-y_l)
    return y
def mollification(c,e,m) :
    # given vectors e and c describing a function and a number m, 
    # mollification two vectors of sized m describing a smoother, 
    # mollified function

    n = e.size 
    if n != c.size :
        print('Vector sizes do not match')

    e_diff = np.zeros(n-1)
    for i in range(0,n-1) :
        e_diff[i] = e[i+1]-e[i]
    epsilon = 3*np.average(e_diff)
    if epsilon < 0 :
        epsilon *= -1

    d = np.arange(-0.9999,0.9999+1.9998/(2*(m-1)),1.9998/(m-1))
    d1 = np.zeros(d.size)
    for i in range(0,d.size) :
        f = - 1/(1-abs(d[i])**2)
        d1[i] = exp(f)
    I_n = 0
    for i in range(0,d.size-1) :
        I_n += (d[i+1]-d[i])*(d1[i]+d1[i+1])/2

    e_min = np.amin(e)
    e_max = np.amax(e)
    e_m = np.arange(e_min,e_max+(e_max-e_min)/(2*(m-1)),(e_max-e_min)/(m-1))
    c_0 = np.zeros(m)
    for i in range(0,m) :
        c_0[i] = piecewise_linear(e_m[i],c,e) 
    c_1 = np.zeros(m) 
    c_m = np.zeros(m)
    for i in range(0,m) :
        x = e_m[i]
        for j in range(0,m) :
            y = e_m[j]
            if abs(x-y) < epsilon :
                
                c_1[j] = c_0[j]*(1/epsilon)*exp(-1/(1-abs((x-y)/epsilon)**2))/I_n
        for j in range(0,m-1) :
            c_m[i] += 0.5*(c_1[j]+c_1[j+1])/(e_m[j+1]-e_m[j])

        c_1 *= 0 
    
    c_min = np.amin(c)
    c_max = np.amax(c)
    c_apparent_height = c_max - c_min
    c_m_min = np.amin(c_m)
    c_m_max = np.amax(c_m)
    c_m_apparent_height = c_m_max - c_m_min
    c_m += - c_m_min
    c_m *= c_apparent_height/c_m_apparent_height
    c_m += c_min
    return c_m,e_m  

def hill_define(peak) :
    # given the coordinates of a peak, hill_define returns 8 points 
    # intended to describe the general area of interest of the hill
    peak_height = dtm.height(peak[0],peak[1])
    t1 = 10100
    t2 = 50
    x = np.arange(peak[0]-t1,peak[0]+t1,t2)
    y = np.arange(peak[1]-t1,peak[1]+t1,t2)
    H = np.zeros([x.shape[0],y.shape[0]])

    for i in range(0,x.shape[0]) :
        for j in range(0,y.shape[0]) :
            H[i,j] = dtm.height(x[i],y[j])  

    h_av = np.average(H)
    h_min = np.amin(H)
    while h_av > peak_height :
        t1 = np.round(0.5*t1)
        t2 = np.round(0.5*t2)
        if t1 < 1 :
            break 
        if t2 < 1 :
            t2 = 1
        x = np.arange(int(peak[0]-t1),int(peak[0]+t1),int(t2))
        y = np.arange(int(peak[1]-t1),int(peak[1]+t1),int(t2))
        H = np.zeros([x.shape[0],y.shape[0]])

        for i in range(0,x.shape[0]) :
            for j in range(0,y.shape[0]) :
                H[i,j] = dtm.height(x[i],y[j])  
        h_av = np.average(H)
        
    if h_av < peak_height :
        app_h = peak_height-h_av
        h1 = h_av + 0.6*app_h
    else :
        app_h = peak_height-h_min
        h1 = h_min + 0.8*app_h

    l_max = app_h*3
    l_step = 5
    j_max = l_max/l_step

    d = np.array([[1,0],[0.7071067811865476,0.7071067811865476],[0,1],\
    [-0.7071067811865476,0.7071067811865476],[-1,0],\
    [-0.7071067811865476,-0.7071067811865476],[0,-1],\
    [0.7071067811865476,-0.7071067811865476]])
    b = np.zeros([8,2])   
    p = np.array([peak[0],peak[1]]) 

    for i in range(0,8) :
        h0 = peak_height 
        j = 0
        while h0 > h1 :
            j += 1
            p = np.array([int(np.round(peak[0] + l_step*j*d[i,0])),\
            int(np.round(peak[1] + l_step*j*d[i,1]))])
            h0 = dtm.height(p[0],p[1])
            if h0 > peak_height :
                break
            if j >= j_max :
                break
        b[i,0] = int(p[0])
        b[i,1] = int(p[1]) 
    return b

def landmark_plot(p,e,c,landmarks_coor,landmark_names,title,x,y,H) :
    # plots c against e with title as a title 
    # if a landmark is visible from p, then it is also plotted with
    # its name annotated
    import matplotlib.pyplot as plt
    k = landmarks_coor.shape[0]
    landmarks_coor = np.reshape(landmarks_coor,(int(k/2),2))
    landmarks_z = np.zeros(int(k/2))
    alpha = np.zeros(int(k/2))
    beta = np.zeros(int(k/2))

    c_max = np.amax(c)
    c_min = np.amin(c)
    apparent_height = c_max - c_min

    plt.plot(e,c)
    plt.title(title)
    plt.xlabel('Horizontal Angle')
    plt.ylabel('Vertical Angle')
    plt.xlim(phi+theta/2,phi-theta/2)
    plt.ylim(c_min,c_min+theta/2)
    plt.fill_between(e,c,c_min - apparent_height,facecolor='green')

    for i in range(0,int(k/2)) :
        landmarks_z[i] = dtm.height(int(landmarks_coor[i,0]),int(landmarks_coor[i,1]))
        q = np.array([landmarks_coor[i,0],landmarks_coor[i,1],landmarks_z[i]+10])
        a = q-p
        alpha[i] = np.arctan2(a[1],a[0])
        beta[i] = np.arctan2(a[2],sqrt(a[0]**2+a[1])**2)
        if intervis(p,q,x,y,H) :
            if e[n-1] <= alpha[i] <= e[0] :
                plt.plot(alpha[i],beta[i],'ro')
                plt.annotate(landmark_names[i],[alpha[i],beta[i]])
            elif e[n-1] <= alpha[i] + 2*np.pi <= e[0] :
                plt.plot(alpha[i]+2*np.pi,beta[i],'ro')
                plt.annotate(landmark_names[i],[alpha[i]+2*np.pi,beta[i]])
            elif e[n-1] <= alpha[i] - 2*np.pi <= e[0] :
                plt.plot(alpha[i]-2*np.pi,beta[i],'ro')
                plt.annotate(landmark_names[i],[alpha[i]-2*np.pi,beta[i]])

    plt.show()