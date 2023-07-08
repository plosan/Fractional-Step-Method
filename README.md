# Fractional Step Method

## Navier-Stokes Equations

The Navier-Stokes equations are a set of partial differential equations (PDEs) that describe the motion of a fluid. When one considers an incompressible and constant viscosity fluid, these read

<p align="center">
  <img src="./mdimages/nseqs.svg">
</p>

the first equation being the conservation of mass, and the second equation being the conservation of momentum. Here, density <img src="./mdimages/density.svg"> and dynamic viscosity <img src="./mdimages/viscosity.svg"> are constants; velocity and pressure are functions of position and time, that is, <img src="./mdimages/velocity_function.svg"> and <img src="./mdimages/pressure_function.svg">.

<p align="center">
  <img src="./mdimages/velocity_pressure_function.svg">
</p>


In general, given some initial and boundary conditions, it is not possible to find a solution

## Fractional Step Method

The Fractional Step Method is a computational procedure to solve the incompressible Navier-Stokes equations. Its simplicity and the fact that in general it runs fast, makes it preferrable over other similar methods. 


## Lid-Driven Cavity


