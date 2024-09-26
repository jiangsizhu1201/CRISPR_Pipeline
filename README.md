# CRISPR-jamboree3
Jamboree demonstrating and testing the functionalities of the single Cell Perturb-Seq Pipeline


# Perturbo


To use Perturbo as your guide inference method, you need to create a separate Conda environment for it. 

Follow these steps to set up the environment:

```
### Activate micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate

### Create a micromamba environment
micromamba create -n perturbo_env

### Activate the perturbo_env
micromamba activate perturbo_env

### Install Perturbo using the wheel file (located in the IGVF_Workflows directory)
pip install perturbo-0.0.1-py3-none-any.whl

### Return to the base environment after installation
micromamba deactivate

```
