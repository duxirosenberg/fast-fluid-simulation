import matplotlib.pyplot as plt
import numpy as np

"""
This file is a translation of the baseline from latticeboltzman.py to 
a c style implementation, which does not use numpy. 

This code is not actually used, but this served as a way to understand
the baseline code and how it works.
""" 

def main():
    """ Lattice Boltzmann Simulation """
    
    "Simulation parameters"
    Nx                     = 400    # resolution x-dir
    Ny                     = 100    # resolution y-dir
    rho0                   = 100    # average density
    tau                    = 0.6    # collision timescale
    Nt                     = 5     #4000 # number of timesteps
    
    # Lattice speeds / weights
    NL = 9
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
    weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
    

    """Initial Conditions"""
    F = [1 for i in range(Ny * Nx * NL)] 
    cylinder = [0 for i in range(Nx*Ny)]

    np.random.seed(42)
    noise_arr = 0.01*np.random.randn(Ny,Nx,NL).flatten()
    for i in range(Nx*Ny*NL):
        F[i] += noise_arr[i]

    # Apply the initial conditions
    for y in range(Ny):
        for x in range(Nx):
            index = y * Nx * NL + x * NL + 3
            F[index] += 2 * (1 + 0.2 * np.cos(2 * np.pi * x / Nx * 4))

    rho = [0 for _ in range(Nx * Ny)]  
    for y in range(Ny):
        for x in range(Nx):
            rho_index = y * Nx + x
            for l in range(NL):
                index = y * Nx * NL + x * NL + l
                rho[rho_index] += F[index]

    # Apply the adjustment for each layer based on rho
    for y in range(Ny):
        for x in range(Nx):
            for l in range(NL):
                index = y * Nx * NL + x * NL + l
                rho_index = y * Nx + x
                F[index] *= rho0 / rho[rho_index]

    """Cylinder Boundary"""
    cylinder = [0 for i in range(Nx*Ny)]
    for i in range(Ny):
        for j in range(Nx):
            if (j - Nx/4)*(j - Nx/4) + (i - Ny/2)*(i - Ny/2) < (Ny/4)*(Ny/4):
                cylinder[j+i*Nx] = 1
    
    
    """Simulation Main Loop"""
    for it in range(Nt):
        print(it)

        """DRIFT"""
        temp_F = F[:]
        # Drift in the x direction for all layers
        for l in range(NL):
            for y in range(Ny):
                for x in range(Nx):
                    source_index = y * Nx * NL + x * NL + l
                    shifted_x = (x + cxs[l]) % Nx
                    target_index = y * Nx * NL + shifted_x * NL + l
                    temp_F[target_index] = F[source_index]

        # Copy the results of the x drift to apply the y drift correctly
        F[:] = temp_F[:]

        # Drift in the y direction for all layers, using the updated F from the x drift
        for l in range(NL):
            for x in range(Nx):
                for y in range(Ny):
                    source_index = y * Nx * NL + x * NL + l
                    shifted_y = (y + cys[l]) % Ny
                    target_index = shifted_y * Nx * NL + x * NL + l
                    temp_F[target_index] = F[source_index]

        # Update F with the final results after both drifts
        F[:] = temp_F[:]


        """Set Reflective Boundaries"""
        rearrange_idx = [0, 5, 6, 7, 8, 1, 2, 3, 4]  # The rearrangement pattern
        bndryF = [] 

        for y in range(Ny):
            for x in range(Nx):
                if cylinder[y * Nx + x]:  # If cylinder at (y, x) is True (non-zero)
                    # Select and rearrange elements corresponding to (y, x) in F
                    base_idx = (y * Nx + x) * 9  # Base index for the (y, x) slice in F_flat
                    rearranged = []
                    for i in rearrange_idx:
                        rearranged.append(F[base_idx + i])
                    bndryF.extend(rearranged)

        """Calculate Fluid Variables"""
        rho = [0] * (Nx * Ny)  
        ux = [0] * (Nx * Ny)   
        uy = [0] * (Nx * Ny)   

        for y in range(Ny):
            for x in range(Nx):
                sum_f = 0  
                sum_fx = 0 
                sum_fy = 0 
                
                for l in range(NL):
                    index = (y * Nx + x) * NL + l 
                    f_val = F[index]  # The F value for the current cell and layer
                    
                    sum_f += f_val
                    sum_fx += f_val * cxs[l]
                    sum_fy += f_val * cys[l]
                
                cell_index = y * Nx + x  
                rho[cell_index] = sum_f
                ux[cell_index] = sum_fx / sum_f  
                uy[cell_index] = sum_fy / sum_f  



        """Apply Collision"""
        Feq = [0] * (Nx * Ny * NL)
        for y in range(Ny):
            for x in range(Nx):
                idx = y * Nx + x  # Index for the current cell
                for l in range(NL):
                    cx = cxs[l]
                    cy = cys[l]
                    w = weights[l]
                    cu = 3 * (cx*ux[idx] + cy * uy[idx])
                    uu = ux[idx]*ux[idx] + uy[idx]*uy[idx]
                    Feq[idx * NL + l] = rho[idx] * w * (1 + cu + 0.5 * (cu*cu) - 1.5 * uu)


        for i in range(len(F)):
            F[i] += -(1.0 / tau) * (F[i] - Feq[i])

        """Apply Boundary"""
        bndry_idx = 0  
        for y in range(Ny):
            for x in range(Nx):
                if cylinder[y * Nx + x]:  
                    for l in range(NL):
                        F_index = (y * Nx * NL) + (x * NL) + l
                        # Update F_flat with the corresponding bndryF value
                        F[F_index] = bndryF[bndry_idx]
                        bndry_idx += 1  # Move to the next element in bndryF for the next iteration



if __name__== "__main__":
  main()

