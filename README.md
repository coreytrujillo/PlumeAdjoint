# Tangent Linear and Adjoint Plume Advection Models

## Overview

This repository contains files which demonstrate the adjoint (AD) and tangent linear (TL) models of a basic plume function.

The `plume3.m` function is a version of the advection plume model I have developed with Dr. Daven Henze over the first year of my Ph.D. In terms of adjoint vernacular, this is the "forward model" (FWM).  I simplified this code from the original `plume.m` function to be easier to derive for the AD and TL models and included a cost function. This function was tested against the `plume2.m` function and showed identical results.  

## Codes

The `plume3TLM.m` function is the tangent linear version of `plume3.m`. This is effectively the derivitave of the FWM (`plume3.m`).

The `plume3ADM.m` function is the adjoint of `plume3.m`. This is effectively the TLM in reverse order - it is best to think about this as computationally reversed order rather than temporal or spatial reverse order. Both temporal and spatial sequences are present in the code so thinking in terms of reversing computational steps rather than temporal or spatial steps simplifies things.

The `plume_adjoint_comparison.m` code uses `plume3.m`, `plume3TLM.m` and `plume3ADM.m` to evaluate the gradient of the cost function with respect to emissions (dJ/dE) finite difference (FD), TL, and AD models.

## Background and Details

The FD approach uses two runs of the FWM at each time step. The difference is taken between a base case and a perturbed case with an input emission vector with a perturbation at a single timestep. This difference is used to calculate the gradient of the cost function at that time step. 

The TL approach uses the TLM with a dE vector with a 1 in the current time step and zeros elsewhere. From this, dJ/dE is calculated for that particular time step.

A loop must be written around the FD and TL models to get the gradient of the cost function at each time step. 

The AD approach uses a single perturbed dE value to get the gradient of the cost function at all times.


