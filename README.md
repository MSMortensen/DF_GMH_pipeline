# Created:  2022-05-20	masmo
# Updated:	2022-05-30	masmo	Modified order of commands to install WSL
##############################################################################
###                          GENERAL INFORMATION                         ###
##############################################################################
To run the pipeline you need use a linux environment with conda installed.
On Windows I recommend using Windows subsystem for linux, but a virtual 
machine running Linux would also be fine.

##############################################################################
###                  INSTALL WINDOWS SUBSYSTEM FOR LINUX                   ###
##############################################################################
 
   - Open PowerShell as Administrator and run these two commands:
       dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
       dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
   - restart your computer

   - Go to Microsoft store
   - Install "Ubuntu 20.04.4 LTS" (a different version should be ok)
   - Open Ubuntu and choose username + password
           NOTE: I recommend using your DTU username, but anything works
   - Run the following command:
       sudo apt update && sudo apt upgrade
           NOTE: You need to use the password chosen above

##############################################################################
###                             INSTALL CONDA                              ###
##############################################################################
 
   - Open Ubuntu (or another version of linux)
   - download installation script
        curl -sL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > "Miniconda3.sh"
   - Install Miniconda by entering:
       bash Miniconda3.sh   
   - Restart your Terminal. Now your prompt should list which environment is active
       NOTE: This means that the command prompt begins with: (base)
   - Update Conda using the command:
       conda update conda
   - Install wget for Conda
       conda install wget

##############################################################################
###                            GROUP PIPELINE                              ###
##############################################################################

 All other files necessary to use the pipeline can be copied from:
"O:\GutMicro\Tarmmikrobiologi gruppen - Gut Ecology group\Methods\NGS pipeline & QIIME\Current pipeline"

Further two step are required to finish the general setup of the pipeline:
   1) Install a conda environment by runnning:
       conda env create -f CONDA_env_setup.yml
   2) Copy the folder "DB" to a location that you access when running conda
       After copying the folder go to in and use this command to confirm it:
           echo $PWD
 
                          Now everything is ready.

# To run the pipeline do as follows:
   - copy files in the folder "Code" into your analysis folder
   - Open Lab_16S_Pipeline_DF-22_1.sh and follow the instructions
