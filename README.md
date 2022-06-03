# GENERAL INFORMATION

To run the pipeline you need use a linux environment with conda installed.
On Windows I recommend using Windows subsystem for linux, but a virtual
machine running Linux would also be fine.

# INSTALL WINDOWS SUBSYSTEM FOR LINUX

- Open PowerShell as Administrator and run these two commands:

```
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

- Restart your computer
- Go to Microsoft store
- Install "Ubuntu 20.04.4 LTS" (a different version should be ok)
- Open Ubuntu and choose username + password

> NOTE: I recommend using the same username as your profile name on your computer, but anything works

- Run the following command:

```
sudo apt update && sudo apt upgrade
```

> NOTE: You need to use the password chosen above

# INSTALL CONDA

- Open Ubuntu (or another version of linux)
- download installation script

```
curl -sL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > "Miniconda3.sh"
```

- Install Miniconda by entering:

```
bash Miniconda3.sh
```

- Restart your Terminal. Now your prompt should list which environment is active

> NOTE: This means that the command prompt begins with: (base)

- Update Conda using the command:

```
conda update conda
```

- Install wget for Conda

```
conda install wget
```

# GROUP PIPELINE

All files necessary to use the pipeline can be downloaded from github.

To do this open Ubuntu and run the following command:

```
git clone https://github.com/MSMortensen/DF_GMH_pipeline
```

Then enter the folder and delete unnecessary files

```
cd DF_GMH_pipeline
rm -r .git .gitignore
```

Install the necessary conda environment by runnning:

```
conda env create -f CONDA_env_setup.yml
```

### To run the pipeline do as follows:

- Create a copy of the **analysis** folder and use that folder for your analysis
```
cp analysis ANALYSIS_FOLDER
```

- Enter the chosen folder, open GMH_16S_Pipeline.sh and follow the instructions
```
cd ANALYSIS_FOLDER
```
