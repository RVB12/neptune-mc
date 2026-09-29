# This is a very 1D simple advection diffusion equation solver using a finite volume method.
# It is used as an extra verification of the Monte Carlo code.

import numpy as np
import scipy.sparse
import scipy.sparse.linalg

def minmod(a, b):
    return np.where(a*b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)

# Calculate the time derivative of the density u for the advection-diffusion equation using a finite-volume method
# The equation is given by:
# du/dt = f(u) = d/dx (-v * u + Dr * d/dx( Df * u )) - R_i * u + S
# Where u is the density, v is the drift velocity,
# Dr is the diffusion rate, Df is the diffusion factor,
# R_i is the ionisation rate, and S is the source term
#
# u0: size n array of initial density in each cell
# v: size n array of drift velocity in each cell
# Dr: size n array of diffusion rate in each cell
# Df: size n array of diffusion factor in each cell
# Ri: size n array of ionisation rate in each cell
# S: size n array of source term in each cell
# xb: size n+1 array of cell boundaries
# bcL, bcR: Coefficiens for the left and right boundary conditions
#           [alpha, beta, gamma]
#           alpha * u_face + beta * du/dx_face = gamma
def advection_diffusion_dudt(u0, v, Dr, Df, Ri, S, xb, bcL, bcR):
    n = len(u0)

    # Cell center of each cell
    # This is simply equal to the boundary values at the boundaries
    xc = np.zeros(n+2)
    xc[1:-1] = 0.5 * (xb[1:] + xb[:-1])
    xc[0] = xb[0]
    xc[-1] = xb[-1]

    # Width of each cell
    dx = np.zeros(n+2)
    dx[1:-1] = xb[1:] - xb[:-1]
    dx[0] = xc[1] - xc[0]
    dx[-1] = xc[-1] - xc[-2]

    # 2 ghost cells are included here for the boundary conditions
    # these ghost cells are the value at the boundary, not some cell outside the domain
    u = np.zeros(n+2)
    u[1:-1] = u0
    # Calculate the ghost cell values based on the boundary conditions
    u[0] = (bcL[2] - (bcL[1] / dx[0]) * u[1])/(bcL[0] - bcL[1] / dx[0])
    u[-1] = (bcR[2] + (bcR[1] / dx[-1]) * u[-2])/(bcR[0] + bcR[1] / dx[-1])

    # Extend the diffusion factor to include ghost cells
    Df_ext = np.zeros(n+2)
    Df_ext[1:-1] = Df
    Df_ext[0] = Df[0]
    Df_ext[-1] = Df[-1]

    # Interpolate coefficients to cell faces
    # The i-th face is between the i-1 and i cell centers
    Dr_face = np.zeros(n+1)
    Dr_face[1:-1] = 0.5 * (Dr[:-1] + Dr[1:])
    Dr_face[0] = Dr[0]
    Dr_face[-1] = Dr[-1]

    v_face = np.zeros(n+1)
    v_face[1:-1] = 0.5 * (v[:-1] + v[1:])
    v_face[0] = v[0]
    v_face[-1] = v[-1]

    # MUSCL reconstruction for advection term with minmod slope limiter
    # slopes at cell centers (size n+2, but only valid for 1..n)
    slope = np.zeros_like(u)
    dxm = xc[1:-1] - xc[:-2]   # x_i - x_{i-1}
    dxp = xc[2:]   - xc[1:-1]  # x_{i+1} - x_i
    slopeL = (u[1:-1] - u[:-2]) / dxm
    slopeR = (u[2:]   - u[1:-1]) / dxp
    slope[1:-1] = minmod(slopeL, slopeR)
    # --- Reconstruction at faces ---
    dxL = xb - xc[:-1]   # distance from left cell center to face
    dxR = xc[1:] - xb    # distance from right cell center to face
    # left state at face i (from cell i)
    uL = u[:-1] + slope[:-1] * dxL
    # right state at face i (from cell i+1)
    uR = u[1:] - slope[1:] * dxR
    # upwind selection
    u_face = np.where(v_face > 0, uL, uR)
    u_face[0] = u[0]  # left boundary face value
    u_face[-1] = u[-1]  # right boundary face value
    # flux
    adv_flux = v_face * u_face
    # adv_flux = np.where(v_face > 0, v_face * u[:-1], v_face * u[1:])

    # Second order central difference for diffusion term
    diff_flux = np.zeros(n+1)
    diff_flux = -Dr_face * (Df_ext[1:] * u[1:] - Df_ext[:-1] * u[:-1]) / (xc[1:] - xc[:-1])

    dudt = (-adv_flux[1:] + adv_flux[:-1] -diff_flux[1:] + diff_flux[:-1]) / dx[1:-1]
    dudt += -Ri * u[1:-1] + S

    return dudt