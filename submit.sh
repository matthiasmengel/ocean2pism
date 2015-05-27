#!/bin/bash

# This file is part of ocean2pism.

# ocean2pism is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ocean2pism is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ocean2pism.  If not, see <http://www.gnu.org/licenses/>.


# Run this script from your shell to submit to the cluster.
# Choose the number of cores with -t.
# The program causes heavy file I/O when used with many cores,
# i'd suggest to use not more than 48.
# Add export PATH=$PATH:/iplex/01/sys/applications/python/ibin to your .profile or .bashrc

#llrun -t 4 -c medium -g isimip -N Diffusion cpy runDiffusion.py
llrun -t 50 -c short -g isimip -N Diffusion cpy runDiffuse.py

