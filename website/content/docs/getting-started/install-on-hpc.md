---
title: "Install and run on HPC"
description: "Install on HPC"
lead: ""
date: 2021-05-26
lastmod: 2021-05-26
draft: false
images: []
menu:
  docs:
    parent: "getting-started"
weight: 40
toc: true
---

This page describes how to install and run Eilmer on
some of the HPC cluster systems our research group
typically has access to.
If you have access to another cluster, you still might
find the information here useful.
You might find your local cluster setup is similar
to one of those listed here.
Likely, you should be able to adapt the instructions with
minimal hassle.

This section assumes you know how to set up your environment
for running eilmer on the cluster machine. What we give
you here are specific steps required to build and run
the MPI version.

> **A note on D compiler versions**
>
> As D is a relatively new language, and particularly so
> in HPC environments, we find it best to bring-your-own
> D compiler. We recommend having a very recent version
> of the D compiler installed in your own account and
> commit to updating the compiler on a regular basis.
> For optimised builds use the LLVM compiler.
> Install notes for the LLVM D compiler (LDC) are
> available [here]({{< relref "install-d-compiler" >}}).


## A 3-minute introduction to queue systems
HPC systems are shared access systems for multiple users.
The compute jobs are controlled by a queueing system.
Different queueing systems are installed on various clusters,
but they tend to have some commonalities.
Those commonalities are: a way to submit your jobs to the queue;
a way to check on your jobs' status; and a way to remove
a job from the queue or terminate it while it's running.

In the clusters we use regularly, there are two queueing
systems: PBS(/Pro) and SLURM.
Here is a brief summary on how to interact with those systems.

### Preparing a submission script
Your compute job will run in batch mode on the cluster.
What this means is that you prepare the instructions (ie. the list of commands you would type)
ahead of time and place these in a text file.
This is what we call the submission script.
We also include some directions to the queue system itself in the submission script.
So, there are two parts in a submission script:

  1. Directives for the queue system. Each of these lines start with `#` symbol
     and a phrase associated with the particular system. On PBS, use `#PBS`.
     On SLURM, use `#SBATCH`.
  2. Commands to launch your compute job. These are the commands as you would
     actually type them if running your job interactively at a terminal.

> We need to add a caveat on that last statement.
> You do type commands as you would at the terminal to start your job,
> but remember your job is running in batch mode and you won't see its
> output or error messages.
> For this reason, it's common to redirect that output with the stdout redirect `>`
> and stderr redirect `2>`.
> You will find this in the examples below.

Here is an example script to give those points above some concreteness.


    #!/bin/bash
    #SBATCH --job-name=my-MPI-job
    #SBATCH --nodes=1
    #SBATCH --ntasks=24

    e4shared --job=myJob --run  > LOGFILE 2> ERRFILE

The first line sets our shell.
The next three lines are directives for the SLURM queue manager.
We won't explain them here because they are specific to each HPC cluster.
Here, we are just emphasising the general layout of a submission script.
There is one command in this script: a command to launch Eilmer.
Note the stdout is directed to a file called `LOGFILE` and
stderr to `ERRFILE`.

Save your text file submission script.
The extension on the script doesn't matter.
I usually save PBS scripts with a `.qsub` extension
and SLURM scripts with a `.sbatch` extension.
The reason is related to submission: those extensions remind we
which queue system I prepared my job for.

### Submitting a job

After saving your submission script into a file,
you are ready to submit this to the queue system.
The submission command in PBS is `qsub`.
In SLURM, use `sbatch` to submit a job.
This typed at the command prompt on the login node of a cluster.
For example, assuming a PBS system and submission script
called `run-my-job.qsub`, type:

    qsub run-my-job.qsub

The system will return to you a job number.

### Checking job status
Your job may not launch instantly.
In fact, it might be several hours before your job starts.
To check the job's status (or multiple jobs, if you have launched multiples),
use the queue status command.
On PBS, use `qstat`.
On SLURM, use `squeue`.

### Removing a job
If your job runs successfully to completion, you don't need to do anything.
The job simply terminates and leaves output in your working directory.

Occasionally, you might detect an error in your script or
you might have changed your mind about parameters in your job.
You will need a way to remove your job from either the queue,
or running on the cluster.
Here is where your job number comes in handy.
You can remove a job using the appropriate command followed
by your job number.
On PBS, delete a job using `qdel`.
On SLURM, cancel a job using `scancel`.

This example shows what to type to remove a job with a
job id `3465` on a SLURM system:

    scancel 3465


> As tempting as it might be to remove others' jobs so your
> job gets on faster, these commands won't let you do that.
> You can only remove your jobs from the queue.


### Summary of commands and where to find more information

There are several more commands available to interact
with a queue system.
However, the main three you will need are: a command to submit a job;
a command to check status; and a command to remove a job if you need.
These are summarised here.

 | Action           | `PBS`              | `SLURM`                |
 |------------------|--------------------|------------------------|
 | submit a job     | `qsub my-job.qsub` | `sbatch my-job.sbatch` |
 | check job status | `qstat`            | `squeue`               |
 | remove a job     | `qdel 3654.pbs`    | `scancel 3654`         |

There is more information that is specific to each of the clusters
in the sections below.
You should consult the user guides for specific clusters for hints
on the queue submission.
You can also use `man` to find out the full list of options
available on these queue commands.
What I've introduced here is just their very basic usage.



## Goliath: EAIT Faculty Cluster

*Hardware:* Dell servers, Intel CPUs, 24-cores per node

*Operating system:* CentOS Linux 7.4

*Queue system:* SLURM

### Compiling
As this is a modern RedHat-flavoured system, you will need
to load the openmpi module to both build and run `e4mpi`.

To compile:

    module load mpi/openmpi-x86_64
    cd dgd/src/eilmer
    make FLAVOUR=fast WITH_MPI=1 install

To complete the install, remember to set your environment variables
as described on the [install page.]({{< relref "installation-locally#setting-environment-variables" >}})


### Running

Where possible, try to use a full node at a time on Goliath.
That means trying to set up your simulations to use 24 cores
per node.
Let's look at an example SLURM submit script where you
have 24 MPI tasks to run on 1 node.
I've named this file `run-on-goliath.slurm`.

    #!/bin/bash
    #SBATCH --job-name=my-MPI-job
    #SBATCH --nodes=1
    #SBATCH --ntasks=24

    module load mpi/openmpi-x86_64
    mpirun e4mpi --job=myJob --run  > LOGFILE

We submit that job using `sbatch`:

    sbatch run-on-goliath.slurm

Next, we'll look at running an oversubscribed job.
An oversubscribed job is when we setup *more* MPI tasks
than we have cores.
This might be the case if we have many blocks (eg. 120) and we
wish to run the 24 cores of a single node.


    #!/bin/bash
    #SBATCH --job-name=my-oversubscribed-MPI-job
    #SBATCH --nodes=1
    #SBATCH --ntasks=120
    #SBATCH --overcommit

    module load mpi/openmpi-x86_64
    mpirun e4mpi --job=myJob --run  > LOGFILE

Note the main changes are to correctly specify the number of
MPI tasks (here we have 120) and to add the directive `--overcommit`.


## Tinaroo: UQ RCC Cluster

*Hardware:* SGI severs, Intel CPUs, 24 cores-per-node, InfiniBand interconnect

*Operating system:* Rocks using CentOS Linux 6.9

*Queue system:* PBS

### Compiling
Tinaroo has an extensive module system with various options for
compilers and MPI libraries.
These have been compiled specifically on tinaroo to take
advantage of its InfiniBand interconnect.
Select an openmpi environment for compiling.
The `openmpi2_ib/2.1.6` has been tested and works nicely with the code.
To compile Eilmer4 with MPI enabled on tinaroo, do:

    module purge
    module load gnu
    module load openmpi2_ib/2.1.6
    cd dgd/install-scripts
    ./install-transient-solvers.sh

To complete the install, remember to set your environment variables
as described on the [install page.]({{< relref "installation-locally#setting-environment-variables" >}})

### Running
Where possible, try to use a full node at once on tinaroo.
The following is a tinaroo submit script that will request
24 MPI tasks. We have set up our job such that it is split
across 24 MPI tasks.

    #!/bin/bash
    #PBS -S /bin/bash
    #PBS -N my-MPI-job
    #PBS -A UQ-EAIT-MechMining
    #PBS -l select=1:ncpus=24:mpiprocs=24:mem=10g
    #PBS -l walltime=24:00:00

    module purge
    module load gnu
    module load openmpi2_ib/2.1.6

    cd $PBS_O_WORKDIR
    mpirun e4mpi --job=dbl-cone --run > LOGFILE

The `select` value is used to select the number of nodes.
The remaining values related to values per node.
So we've requested 24 CPUs on a node and 10 GB of memory
on a node.
If you set `select=3`, then the total resource request
would be for 72 CPUs and 30 GB of memory.
We use `qsub` to submit our job script:

    qsub run-on-tinaroo.qsub

## Gadi: NCI Cluster

*Hardware:* Fujitsu servers, Intel 'Cascade Lake' processors: 48 cores-per-node (two 24-core CPUs), 192 Gigabytes of RAM per node, HDR InfiniBand interconnect

*Operating system:* Linux, CentOS 8

*Queue system:* PBSPro

### Alternative to Compiling: Using the pre-built module

Gadi has a number of helpful modules already present in the /scratch/dc7/common directory, including a copy of Eilmer 4.
To avoid having to compile the code yourself, it is possible to simply load this module and allow it to configure your environment variables.
First set up  privatemodules directory in your home directory.

    $ cd
    $ mkdir privatemodules
    $ cd privatemodules
    $ mkdir dgd
    $ ln -s /scratch/dc7/common/dgd/ldc2/dgd-module dgd/ldc

Then add the following lines to the .bash\_profile file in your home directory.
Note that this file is preferred over modifying .bashrc on Gadi, since it is always loaded when you start a new shell.

    module load gcc
    module load openmpi/2.1.6
    module load dgd/ldc

Finally, log out of Gadi and back in again to load the new modules. You can check the paths have been set properly by typing:

    $ which e4shared
    /scratch/dc7/common/dgd/ldc2/bin/e4shared

    $ e4shared version
    Eilmer 4.0 compressible-flow simulation code.
    Revision: 250f175e
    Compiler-name: ldc2
    Build-flavour: fast
    Capabilities: multi-species-gas multi-temperature-gas MHD turbulence.
    Parallelism: Shared-memory

### Compiling

If you require your own copy of the source code, you will need to compile Eilmer yourself.
Since Gadi does not have a module for the ldc2 compiler, you will need to install it too: See the instructions in "Guide to installing LDC".
(There is also a pre-built ldc compiler in /scratch/dc7/common.)
In terms of modules, you only need to load the openmpi module.
I have had success with the 4.0.2 version.

    module load openmpi/4.0.2

As described for other systems, use the `install-transient-solvers.sh` script to
get an optimised build of the distributed-memory (MPI) transient solver.

    cd dgd/install-scripts
    ./install-transient-solvers.sh

To complete the install, remember to set your environment variables
as described on the [install page.]({{< relref "installation-locally#setting-environment-variables" >}})

### An example setup in `.bash_profile`
Here, I include what I have added to the end of my `.bash_profile` file on Gadi.
It sets my environment up to compile and run Eilmer.
Note also that I configure access to a locally installed version of the ldc2 compiler.

    export DGD=${HOME}/dgdinst
    append_path PATH ${DGD}/bin
    append_path PATH ${HOME}/opt/ldc2/bin
    append_path PYTHONPATH ${DGD}/lib

    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so

    module load openmpi/4.0.2

### Running

As is common on large cluster computers, you will need to request the entire node CPU resources if your job spans multiple nodes. On Gadi, that means CPU request numbers are in multiples of 48.
Here is a submit script where I've set up 192 MPI tasks for my job.

    #!/bin/bash
    #PBS -N my-MPI-job
    #PBS -P dc7
    #PBS -l walltime=00:30:00
    #PBS -l ncpus=192
    #PBS -l mem=200GB
    #PBS -l wd
    #PBS -l storage=scratch/dc7

    mpirun e4mpi --job=myJob --run > LOGFILE

Jobs on Gadi are submitted using `qsub`.

Take note about the `storage` directive that appears in the submission script.
It is strongly encouraged on Gadi to set up your jobs in the `scratch` area associated with your
project.
On Gadi, you need to request explicit access to the `scratch` area in your submission
script so that filesystem is available to your job.
The `storage` directive used to request that access.

In the next example, I have set my job up to run in an oversubscribed mode.
Here I have 80 MPI tasks but I'd like to use only 48 CPUs.

    #!/bin/bash
    #PBS -N my-oversubscribed-MPI-job
    #PBS -P dc7
    #PBS -l walltime=00:30:00
    #PBS -l ncpus=48
    #PBS -l mem=50GB
    #PBS -l wd
    #PBS -l storage=scratch/dc7

    mpirun --oversubscribe -n 80 e4mpi --job=myJob --run > LOGFILE-oversubscribed

