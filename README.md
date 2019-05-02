# hyper-planes-in-images
Minimum density hyper planes utilised for image segmentation.

This application allows one to interactly segment images using hyperplanes.  One can either manually select a decision rule by selecting a position on the density or automatically using projection pursuit in accordance of obtaining a minimum density hyperplane solution.

Once a decision boundary is selected, we can improve upon the hyerplane solution by reassigning values within a region around the hyperplane (gamma region) using gradient ascent and thus allowing for a final non-linear solution.

This program is a desktop application that runs using RStudio and Shiny.  
To begin, create a new project in RStudio and place the server.R, ui.R and www folder inside said project folder.
Once the files are located within the correct directory, simply open both R files within RStudio project and run the application.
