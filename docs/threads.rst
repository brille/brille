=================
Parallel threads
=================

``brille`` runs its parallel sections on a pool of native threads.
Methods that take a ``threads`` argument use that many threads when it is one
or more. When it is less than one, which is the default, the number comes from
the ``BRILLE_NUM_THREADS`` environment variable, or is one thread per logical
core if that is not set.

.. code-block:: bash

  BRILLE_NUM_THREADS=4 python my_script.py

Set it when several processes use ``brille`` at once, for example parallel fits
or ``pytest-xdist`` workers, so that together they don't ask for more threads
than the machine has cores. The value must be a positive integer; anything else
is ignored.

It is read each time the pool is sized, so it can also be changed from Python
between calls:

.. code-block:: python

  import os
  os.environ["BRILLE_NUM_THREADS"] = "2"

``OMP_NUM_THREADS`` has no effect, since ``brille`` no longer uses OpenMP.
