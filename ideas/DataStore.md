# Data storage #

It would be good to redesing the data storage interface in the grid cell
network so that it is file format independent.

The proposed workflow:
    1. User submits only *one* script in which they define the parameters and
       simulation jobs.

    2. One simulation comprises one or more simulation runs, each simulation
       run will produce one *set of data*. The data set will be transparent to
       the user, i.e. the user only stores identifiers of their data (and
       potentially can have a file-system-like hierarchy, will be decided
       later)

    3. The simulation has a root directory, defined by the user, that contains
       data from all the simulation runs. The root directory is also a unique
       identifier for the simulation. Ideally, the simulation root directory
       contains in its name a date and time, and a short identifier, which the
       user can set up in the script he/she submits.

    4. The user/submitter passes all the necessary information to the actual
       script that performs the simulation (python in our case)

    5. The script instantiates the data storage object and simply passes the
       objects to store, together with their identifiers.

    6. Simulation runs and stores the data into the appropriate *simulation
       root directory*.

    7. The user performs data analysis, by using the same set of objects as in
       point 5. The user passes the root directory, and simulation identifier
       (number) to the constructor of the data storage object, and specifies
       which data to access by the variable identifiers. The data storage
       object returns (a copy of) the appropriate data on disk/other medium.


## How generic should the data store class be? ##

An important thing to decide is how the user sees the hierarchy of one
simulation run (not the whole simulation which can contain one or more
simulation runs).

A simulation run consists of a place holder for the saved objects. How the
place holder will be implemented can depend on these factors:
    * Is the data targeted for a filesystem storage or database storage?
    * Is the hiearchy in one place holder flat or does it need to be
      multi-level?
    * Can we safely assume that the user will rather produce *big* data that
      are more suited for filesystem storage, rather than data more suitable
      for database storage?

One needs to satisfactorily answer these questions before any implementation
commences, otherwise we will run into troubles.


## Data set realization ##
