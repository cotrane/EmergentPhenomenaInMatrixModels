
The different folders contain the code written for the simulation of different matrix models included in the thesis "Emergent Phenomena in Matrix Models", arXiv: .... A copy of the thesis is also included in the folder.

The different .txt and .py files contain the parameters that can be chosen for the various models. The name of the file corresponds to the matrix model that is simulated.


----- linking to libraries --------

The simulation should compile and run when left in the current directory, if you supply the correct paths to the libraries. In order to achieve this execute the following commands (this commands work for Ubuntu! you need to have superuser rights to do this!):


sudo bash -c "echo 'export WRKSPACEdir=path/to/working/directory' >> /etc/profile"  # this adds the path for the working directory to your environment variable; for example "/home/EmergentPhenonmenaInMatrixModels"
sudo bash -c "echo 'export CLAPACKdir=/path/to/directory' >> /etc/profile" # this adds the path to CLAPACK library to the environment file


To update the env variable immediately in this shell you need to execute "source /etc/profile". For the changes to affect any applications launched outside this shell you need to restart your system.

Now you can check if it was succesful by executing "echo $WKRSPACEdir". This should print the path to the workspace directory on the screen.


For the parallelized versions, HMC_MPI and HMCYM_MPI, which use OpenMPI, this library needs to be installed. The version used is OpenMPI 1.4.3 and needs to be installed on the system in order to use this programs


If you do not have sudo rights you can edit the '~/.bashrc' file. BUT the environment variables in this file will only be loaded in the terminal you log into! 

sudo bash -c "echo 'export WRKSPACEdir=path/to/working/directory' >> ~/.bashrc"
sudo bash -c "echo 'export CLAPACKdir=/path/to/directory' >> ~/.bashrc"


----- building order ---------

In order to run the different models you first need to build all the libraries, which are:

RandomGens, MatrixMan, RepsSUn, GammaMatrices, eigenvalues


------ starting the simulation ---------


The various simulation can be started by the command "./NameOfSimulation"; the parameters will then be loaded from the 'NameOfSimulation.txt' file. If more simulations shall be started after each other with varying parameters execute the simulation using a python file and the command "python NameOfSimulation.py". Note that such a file does not exist for all simulations! It could be written in the same way as for the simulations where such a file is present. The necessary code is included in the projects.

For the parallelized simulations HMC_MPI and HMCYM_MPI the command to execute the simulation is, e.g., "mpirun -n 2 ./HMC_MPI". Here, the executable HMC_MPI would be started using 2 cores. This number can be increased up to the maximal number of matrices, 3 for the 3MM and 8 for the 8MM. More cores will not speed up to simulation but might slow it down.



