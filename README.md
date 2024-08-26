## Terminus ablation calculation and analysis
##still in progress!!
For queries, email: [amankc@u.boisestate.edu](mailto:amankc@u.boisestate.edu)

The datasets and methodologies used to construct the terminus ablation time series are illustrated in Fig. 1.\

![Alt text](Pictures/Updated_Main_flowchart.png)
Figure 1. Flowchart describing the methodology used to calculate terminus ablation. Rounded rectangles show inputs, intermediate outputs,
and final output. Angular gray rectangles show the mathematical processes and purple rectangles show the processes applied to satellite
images. Orange rectangles show data from external sources

![Alt text](Pictures/Frontal_Process.png)
Figure 2: Schematic sketch showing glacier terminus ablation. Schematic sketch showing glacier terminus ablation. The green terminus
position shows glacier advance, blue shows glacier retreat and solid black shows no change in terminus. Solid orange line shows a stationary
gate used to calculate ice discharge. Here, ∆M/∆t is change in terminus over time, A_terminus is Terminus Ablation and U is the velocity.

# Running the code

MATLAB files can be used to generate monthly terminus ablation time series and python notebooks are mostly for analysis of the time series and regional comparison. 

Note about the time series: If there are multiple terminus traces in a single month (mostly during summer months), the terminus ablation is weighted. However, in case of no data, linear interpolation is used to fill the data gap based on two near-month observations.

Matlab sections are numbered. Some of the sections are marked optional. You can run that, depending on the data you are importing or if you want to plot figures.

1. Import all the terminus traces, and add any additional traces you have.
2. Filter the trace to avoid any duplication.
3. Read all the DEMs and landsat images.
4. Read the fjord boundaries. If not, it will prompt you to draw one based on the Landsat image you imported on step 3.
5. Get closest point on Fjordwall to end/startpoint of Terminus Position
6. 
