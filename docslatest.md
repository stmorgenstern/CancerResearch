# Cancer Model Research Files

## Overview
On this page, we will review the structure and functionality of our model. We will provide a brief description of each file and the important functions.

## modelanalysis.jl

This file is the primary interface for the research model. In this file you can specify the number of nodes for discretization and the known physiological parameters and it will generate a plot of the pressure vs radius, velocity vs radius, concentration vs radius, and accumulation vs time.

## soluteperm.jl
Given the hydraulic conductivity and the drug particle radius this file will be used to calculate the vascular permeability and the solute reflection coefficient.

## isolatedvelocity.jl and isolatedpressure.jl
Given the spatial discretization and physiological parameters each of these files will numerically integrate an ODE corresponding to the velocity and pressure as a function of the tumor radius.

## isolatedmodel.jl
This file provides the functionality for integrating a spatial and time dependent PDE for concentration of drug particles within the tumor. Within this file, the problem is formulated for two domains: a low peclet regime and a high peclet regime. The low peclet regime assumes mass transport is taking place via diffusion and convection. The high peclet regime assumes that mass transport is only taking place via convection. The underlying differential equation changes for each regime so these models are implemented independently of eachother. 

Specifically, the *Isolated_Model* and *Isolated_Model_HP* functions utilize the DifferentialEquations.jl to solve the PDE for drug concentration over time and space for the low and high peclet regimes respectively. The code for the discretized PDE's can be found within the *MST!* and *MSTHighPeclet!* functions.

The *ExplicitEulerMSTHighPeclet!* function solves the PDE for the high peclet regime using an explicit Euler integration. This function calls *MSTHighPecletEE!* for the discretized PDE model.

The *Accumulation_Model* function takes the solution matrix from one of the previous functions and averages it over the spatial domain. This provides the accumulation of drug particles in the entirety of the tumor over time. If the solution has less than 100 time nodes it will preserve the original dimensionality but if it is over 100 time nodes it will condense it down to 100 nodes.

## makeplots.jl
This file provides some useful functions to handle the solution matrix of the PDE and provide plots using the Plots.jl package.

The *singleconcplot* and *doubleconcplot* functions will return a plot of concentration vs radius at 10 min, 30 min, and 60 min for one and two solutions respectively.

The *errorconcplot* will provide a plot of percent error vs radius at the 10 min, 30 min, and 60 min time nodes for two different concentration solution matrices.

The *singleaccumplot* and *doubleaccumplot* functions will return a plot of accumulation vs time for one and two solutions respectively.

The *Accumerrorplot* will provide a plot of percent error vs time for two different concentration solution matrices.

The *velocityplot* and *pressureplot* will provide a plot of velocity vs radius and pressure vs radius.

## optimization.jl
This module is a parameter estimation problem that solves for effective permeability by minimizing the sum of squared error between an empirical model based on effective permeability and the integrated PDE model. The model is currently formulated in JuMP using the Ipopt optimizer.
