MPP decomposition for land suppression
++++++++++++++++++++++++++++++++++++++

Before doing long runs it is important to optimise MPP decompositoin by invoking
 land supression to save redundant ocean processors.

In the namelist there is a subsection nammpp. This allows the user to switch on
land suppression from the NEMO simulations. Setting jpni, jpnj and jpnij to values
less than 1 tells NEMO to decompose the domain into jpni x jpnj processor maps
such that jpni x jpnj = total number of requested processors. The user can take
a little control by specifying how the x-axis and y-axis are subdivided
e.g. by setting jpni=10, jpnj=10 and jpnij=100 the x and y axes are divided
up into 10 equal parts for the decomposition. If 40% of the domain happens
to be land land suppression can be triggered if jpni x jpni > jpnij. As a
first guess we may run with the following jpni=10, jpnj=10 and jpnij=60. To test
this run the model for 1time step. Invariably this will fail, but the mode will
write out what jpnij it was expecting in the ocean.output file. jpnij can then
be updated in the namelist_cfg and the model rerun with the correct jpnij value.
This allows the user to run a more efficient setup and speeding up simulations.

When choosing jpni and jpnj bear in mind the following:
* Ideally you want to keep the decomposition as ‘square’ as possible to minimise communications
* With NEMO the efficiency drops of if the tile/processor map size is less than 15x15 points

Initially run with the namelist_cfg parameters::

  !-----------------------------------------------------------------------
  &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
  !-----------------------------------------------------------------------
     cn_mpi_send =  'I'      !  mpi send/recieve type   ='S', 'B', or 'I' for standard send,
                             !  buffer blocking send or immediate non-blocking sends, resp.
     nn_buffer   =   0       !  size in bytes of exported buffer ('B' case), 0 no exportation
     ln_nnogather=  .false.  !  activate code to avoid mpi_allgather use at the northfold
     jpni        =  -1       !  jpni   number of processors following i (set automatically if < 1)
     jpnj        =  -1    !  jpnj   number of processors following j (set automatically if < 1)
     jpnij       =  -1    !  jpnij  number of local domains (set automatically if < 1)


Then inspect ``ocean_output`` to see the ``jpni`` and ``jpnj``.
 Set ``jpnij < jpni*jpnj`` and run again. E.g.::

   !-----------------------------------------------------------------------
   &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
   !-----------------------------------------------------------------------
      cn_mpi_send =  'I'      !  mpi send/recieve type   ='S', 'B', or 'I' for standard send,
                              !  buffer blocking send or immediate non-blocking sends, resp.
      nn_buffer   =   0       !  size in bytes of exported buffer ('B' case), 0 no exportation
      ln_nnogather=  .false.  !  activate code to avoid mpi_allgather use at the northfold
      jpni        =  12       !  jpni   number of processors following i (set automatically if < 1)
      jpnj        =  8    !  jpnj   number of processors following j (set automatically if < 1)
      jpnij       =  1    !  jpnij  number of local domains (set automatically if < 1)

Inspect ``ocean_output`` to find ``jpnij``. In my simulation ``jpni=8, jpnj=12 --> jpnij = 93``
Update OCEANCORES in runscript (make sure the ``aprun`` statement is as expected too)::

  vi runscript
  ...
  OCEANCORES=93 # formerly 96

And submit again.
