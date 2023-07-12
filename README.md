# Fractional Step Method

## Navier-Stokes Equations

The Navier-Stokes equations are a set of partial differential equations (PDEs) that describe the motion of a fluid. When one considers an incompressible and constant viscosity fluid, these read

<p align="center">
  <img src="./md_images/nseqs.svg">
</p>

the first equation being the conservation of mass, and the second equation being the conservation of momentum. Here, density <img src="./md_images/density.svg"> and dynamic viscosity <img src="./md_images/viscosity.svg"> are constants; velocity and pressure are functions of position and time, that is,

<p align="center">
  <img src="./md_images/velocity_pressure_function.svg">
</p>

Let <img src="./md_images/omega.svg"> be an open set with smooth boundary <img src="./md_images/boundary.svg"> (or at least <img src="./md_images/class_c1.svg">). Given some initial conditions on velocity and pressure

<p align="center">
  <img src="./md_images/in_cond.svg">
</p>

along with some boundary conditions, one wants to show that there exists a unique solution (in some sense) to the following problem:

<p align="center">
  <img src="./md_images/ivp.svg">.
</p>

A problem of this kind 
In general, this problem is quite hard. However, for engineering purporses, it suffices to obtain an approximate solution via a numerical method.



## Fractional Step Method

The Fractional Step Method is a computational procedure to solve the incompressible Navier-Stokes equations. Its simplicity and the fact that in general it runs fast, makes it preferrable over other similar methods. 


## Application: Driven Cavity

As an application of the Fractional Step Method, I solved the Lid-Driven Cavity problem, which is shown below:


<p align="center">
  <img src="./md_images/driven_cavity.svg">.
</p>

It consists of a rectangular cavity filled with fluid and open at the top. On the left, right and lower walls, the no-slip condition is imposed. On the upper wall there is a fluid flow that drives the fluid whithin the cavity. The behavior of the fluid is modeled via the Navier-Stokes equations. 

For simplicity, <img src="./md_images/l_equals_h.svg">; the physical data is the following: 


<p align="center">
  <img src="./md_images/physical_data.svg">
</p>

thus the Reynolds number is <img src="./md_images/reynolds.svg">. As for numerical data:

<ul>
  <li>The domain is discretized in a <img src="./md_images/129.svg"> node-centered mesh</li>
  <li>Pressures are computed at nodal locations</li>
  <li>Velocities are computed at control volume faces using staggered mesh (to avoid a checker board problem)</li>
  <li>QUICK scheme is used to compute velocities. It works particularly well in low Reynolds number situations</li>
  <li>Gauss-Seidel method is used to solve the linear system</li>
</ul> 

Results are shown below for times $

The simulation runs until <img src="./md_images/t_206.svg">

The plots below show the results: left column - velocity field, right column - pressure





<p align="center">
  <img src="./md_images/plots/vel_005.svg" width="325" />
  <img src="./md_images/plots/pres_005.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_030.svg" width="325" />
  <img src="./md_images/plots/pres_030.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_060.svg" width="325" />
  <img src="./md_images/plots/pres_060.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_090.svg" width="325" />
  <img src="./md_images/plots/pres_090.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_120.svg" width="325" />
  <img src="./md_images/plots/pres_120.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_150.svg" width="325" />
  <img src="./md_images/plots/pres_150.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_180.svg" width="325" />
  <img src="./md_images/plots/pres_180.svg" width="325" />
</P>


<p align="center">
  <img src="./md_images/plots/vel_206.svg" width="325" />
  <img src="./md_images/plots/pres_206.svg" width="325" />
</P>



