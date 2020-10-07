# EKI_geophysics_2020
Example MATLAB scripts for running ensemble Kalman inversion for electrical resistivity tomography (ERT) data to accompany paper (submitted).

*You must cite:* at the moment, please contact the author. 

# Ensemble Kalman inversion (EKI) for geophysics (2D and 3D)

This repo documents script files to run Ensemble Kalman inversion with level set parameterization. Its main advantage is to delinate resistvity fields as zones (under uncertainty) at a modest computational cost. There is an option to allow fixing the zonal resistivity values, or allow heterogeniety within each zone.

[R2] is used as the forward solver of the ERT problem.

### Requirements
- MATLAB(R) 
- MATLAB(R) statistical and machine learning toolbox and parallel computing toolbox
- We have tested on several workstations running Linux and Mac OS. Sould also work on Windows machine but some Unix-specific commands needs revising (in progress). Also remove "wine" before "R2.exe" when running on Windows.
- 2D problems:
  - For Mac OS/Linux machine: You will need to install [wine]. See more specific instructions below.
- 3D problems: ERT forward modelling is run by [E4D] to allow . Make sure you install the [E4D development version]. Unfortunately it only runs for Linux systems. You may try using [cygwin] to run it on Windows systems.

### General steps
1. Create forward problem in [R2] by defining a true resistivity model, electrode configuration, and survey design. Make sure it runs without returning erros. Do not allow polyline clipping in `R2.in`. You may set that in `polyline.txt`. Copy the [R2] input files in the working directory.
1.1 For field problems, make sure you replace `protocol.dat` with field data and comment out the line in `EKI.m`.
2. Set up the EKI problem by changing `EKI.m`and `Set_prior.m` (see details below)
3. Run `EKI.m` until the line *`Inversion()`*. This will take several minutes to several hours.
4. To process and plot the results, run the remainder of `EKI.m` as appropriate, which calls helper functions `vtk_read()` and `plot_vtk_2D()`.

### `EKI.m`
- **tuning** : 
- **tuning** : 
- **tuning** : 
- 

The EKI script here works on a quadrilateral grid that interpolates the grid used in [R2]

### `Set_prior.m`
- **tuning** : 
- **tuning** : 
- **tuning** : 

### Outputs:


### Examples:
 Name of folder | description | figure in paper
 ---|---|---
`acronym` | Shorthand for datasets | Fig. 4

### Inversion notes:
- prior length scales are 7:1 for x:y by default

### Plotting and analyzing results (2D):


`vtk_read()` allows you to read a vtk file generate from R2 (not all vtk formats will work).

### Plotting and analyzing results (3D):
The same idea applies for plotting 3D res

See this help document () to get started.

### FAQ:

| Question | Answer |
| ------ | ------ |
Why are two grids specified in `EKI.m`? | *The inversion is conducted in a rectangular grid, while it evaluates its update with the forward models which run in R2. At the end of each iteration, the estiamted resistivity field is updated from the rectangular grid to the R2 grid. Therefore, your resultant grid is in R2 format.*
What does polyline cropping do and why is it used? | *Learn more about polyline in the [R2] manual*
How big should my EKI be? | *It should be a rectangular bounding box that covers the entire region of interest. In oder words, it should cover the entire region that is cropped by the polyline in the above question.*
How to run EKi in parallel? | That's why we implment it MATLAB. With the parallel toolbox installed, it should do so automatically. The code will still run without it but without parallelism. The MATLAB `parfor` function means the content in the for loop is evaluated not necessarily in order.


### Todos

 - Feel free to suggest to opening an issue

### Relevant repos
 - [ResIPy]: great API and graphical user interface to [R2]. You may find it helpful to generate input files for [R2]
 - [E4D tools]: you may find it helpful
 - []
 

License
----

GPL 3.0 License

**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [R2]: <https://github.com/joemccann/dillinger>
   [E4D]:
   [E4D development version]:
   [ResIPy]:
   [E4D Tools]:
   [wine]:
   [cygwin]:
   
   

