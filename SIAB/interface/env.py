"""this module is for converting common Liunx command to the version
compatible with High Performance Computing (HPC) system or say for
Supercomputer."""
import os
import subprocess
import sys
import json

##############################################
#                                            #
##############################################
def submit(folder, 
           module_load_command,
           mpi_command,
           program_command,
           env = 'local',
           runtype = 'dry') -> str:
    """general submit function with compatibility of HPC systems
    
    Parameters
    ----------
    folder : str
        folder to run the command
    module_load_command : str
        module load command
    mpi_command : str
        mpi command
    program_command : str
        program command
    env : str
        environment to run the command
    runtype : str
        type of running the command, can be 'dry-run' or 'run'
    
    Returns
    -------
    str
        job template generator
    """
    assert runtype in ['dry-run', 'run'], "runtype can only be 'dry-run' or 'run'"
    assert env in ['local', 'hpc'], "env can only be 'local' or 'hpc'"

    jtg = "%s\n"%module_load_command
    jtg += "echo \"present directory: \" `pwd`;\n"
    jtg += "export OMP_NUM_THREADS=1\n"
    jtg += "echo \"OMP_NUM_THREADS:\" $OMP_NUM_THREADS\n"
    jtg += "folder=%s\n"%folder
    jtg += "program_command='%s'\n"%(program_command)
    jtg += "mpi_command='%s'\n"%(mpi_command)
    jtg += "echo \"run with command: $mpi_command $program_command\"\n"
    jtg += "stdbuf -oL $mpi_command $program_command"

    # here = os.path.abspath(os.path.dirname(__file__))
    os.chdir(folder)

    hpc_settings = {"shell": True, "text": True, "timeout": 72000}
    
    # run command backup, recording the command to be run, can be easily restart
    with open('command.json', 'w') as f:
        json.dump({'command': jtg, 'env': env, 'hpc_settings': hpc_settings}, f)
    print(f"running command has been saved to file {os.path.join(folder, 'command.json')}",
          flush=True)
    
    # run the command only if runtype is 'run'
    if runtype == 'run':
        run(command=jtg, env=env, hpc_settings=hpc_settings)

    # os.chdir(here) # back to the original directory
    os.chdir("..")
    return jtg

##############################################
#        basic wrapped linux commands        #
##############################################

def run(command,
        env = "local",
        additional_args = None,
        hpc_settings = None):
    """run command in different environment
    
    Parameters
    ----------
    command : str
        command to be executed
    env : str
        environment to run the command
    additional_args : list
        additional arguments to be added to the command
    hpc_settings : dict
        settings for HPC system
    
    Returns
    -------
    int
        return value of the command
    """
    if additional_args is not None:
        command = " ".join([command, *additional_args])
    if hpc_settings is None:
        hpc_settings = {
            "stdin": subprocess.DEVNULL,
            "shell": True,
            "text": True,
            "timeout": 60
        }

    if env == "local":
        value = os.system(command)
    elif env == "hpc":
        sys.stdout.flush()
        value = subprocess.run(command, **hpc_settings)
        sys.stdout.flush()
    return value

def _rerun_folder(folder):
    """rerun the command in the folder
    
    Parameters
    ----------
    folder : str
        folder to rerun the command
    """
    with open(os.path.join(folder, 'command.json'), 'r') as f:
        command = json.load(f)
    return run(command=command['command'], 
               env=command['env'], 
               hpc_settings=command['hpc_settings'])

def op(operation: str,
       src: str,
       dst: str = "",
       env: str = "local",
       additional_args: list = None,
       hpc_settings: dict = None):
    """operation on files
    
    Parameters
    ----------
    operation : str
        operation to be executed
    src : str
        source file
    dst : str
        destination file
    env : str
        environment to run the command
    additional_args : list
        additional arguments to be added to the command
    hpc_settings : dict
        settings for HPC system
    
    Returns
    -------
    int
        return value of the command
    """
    if additional_args is None:
        additional_args = []
    additional_args = [operation, *additional_args]
    return run(command=" ".join([*additional_args, src, dst]).replace("  ", " "),
               env=env,
               hpc_settings=hpc_settings)

