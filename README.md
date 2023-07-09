# Fractional Step Method

## Navier-Stokes Equations

The Navier-Stokes equations are a set of partial differential equations (PDEs) that describe the motion of a fluid. When one considers an incompressible and constant viscosity fluid, these read

<p align="center">
  <img src="./mdimages/nseqs.svg">
</p>

the first equation being the conservation of mass, and the second equation being the conservation of momentum. Here, density <img src="./mdimages/density.svg"> and dynamic viscosity <img src="./mdimages/viscosity.svg"> are constants; velocity and pressure are functions of position and time, that is,

<p align="center">
  <img src="./mdimages/velocity_pressure_function.svg">
</p>

Let <img src="./mdimages/omega.svg"> be an open set with smooth boundary <img src="./mdimages/boundary.svg"> (or at least <img src="./mdimages/class_c1.svg">). Given some initial conditions on velocity and pressure

<p align="center">
  <img src="./mdimages/in_cond.svg">
</p>

along with some boundary conditions, one wants to show that there exists a unique solution (in some sense) to the following problem:

<p align="center">
  <img src="./mdimages/ivp.svg">.
</p>

In general, this problem is quite hard. However, for engineering purporses, it suffices to obtain an approximate solution via a numerical method.



## Fractional Step Method

The Fractional Step Method is a computational procedure to solve the incompressible Navier-Stokes equations. Its simplicity and the fact that in general it runs fast, makes it preferrable over other similar methods. 


## Application: Driven Cavity

As an application of the Fractional Step Method, I solved the Lid-Driven Cavity problem, which is depicted below:


<p align="center">
  <img src="./mdimages/driven_cavity.svg">.
</p>

It consists of a rectangular cavity filled with fluid and open at the top. On the left, right and lower walls, the no-slip condition is imposed. On the upper wall there is a fluid flow that drives the fluid whithin the cavity. 





