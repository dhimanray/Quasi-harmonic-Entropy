#--------------------------------------------------------------------------------
#This code calcultes the Quasi-harmonic Entropy according to the method proposed
#in the following publication
#
#
#Andricioaei, Ioan, and Martin Karplus. "On the calculation of entropy from 
#covariance matrices of the atomic fluctuations." The Journal of Chemical 
#Physics 115.14 (2001): 6289-6292.
#-------------------------------------------------------------------------------

import numpy as np

import MDAnalysis as md

import MDAnalysis.analysis.encore.covariance as cov

t = np.arange(5000,45000,1000) #choose according to the name of trajevtory files

for time in t:

    u = md.Universe('TATA.pdb','TATA_traj_%d_ps.dcd'%time)
    heavy = u.select_atoms('all')
    frames = len(u.trajectory)

#---------- get the coordinates in a (time x 3N) numpy array --------------------#    
    n_atoms = len(heavy.positions)
    x = np.zeros((frames,n_atoms*3))
    #print(x.shape)
    for ts in u.trajectory:
        y = heavy.positions
        x[ts.frame] = y.flatten()
    #print(len(x[:,1]))

#------------- Subtract Mean ---------------------------------------------------#
    for i in range(n_atoms*3):
        x[:,i] -= np.mean(x[:,i])
    #print(np.mean(x[:,23]))

#---------- Calculate covariance matrix --------------------------------------- #
    c = (1./float(frames))*x.T.dot(x)
    #print(c[1,3],c[3,1])

#---------- Do mass weighting ------------------------------------------------- #
    m = u.atoms.masses
    #print(m)
    m = np.repeat(m, 3)
    m_half = np.diag(m**0.5)
    #print(m)
    sigma = np.dot(m_half,(np.dot(c,m_half)))

#---------- Calculate eigenvalues --------------------------------------------- #
    lam, v = np.linalg.eigh(sigma)
    lam = lam[::-1]

    #print(lam[:])


    #temperature
    T = 298.0
    kB = 0.6/300.0 #kcal/mol
    #calculate entropy
    S = 0.0
    
    f1=open('enmode_%d.dat'%time,'w')
    for i in range(len(lam)):
        wi = np.sqrt(1.0/lam[i])
        #---------------------------------------
        #In SI unit wi = (np.sqrt(1.0/lam[i]))*1.569*1E+13
        #hbar/kBT = 0.2566*1E-13
        #hbar*wi/kBT = 0.4025 * wi
        #----------------------------------------------
        a = 0.4025 * wi   

        #Andricioaei-Karplus formula (JCP 2001)
        S_i = a/(np.exp(a)-1.0) - np.log(1-np.exp(-a))

        #Schlitter formula
        #S_i = 0.5*np.log(1.0 + (np.exp(1.0)/a)**2)

        S += S_i*kB
        print(i,S_i*kB,file=f1)
    
    print(float(time), S)


