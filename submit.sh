
# Run this script from your shell to submit your IMOA calc to the cluster.
# Choose the number of cores with -t.
# The program causes heavy file I/O when used with many cores,
# i'd suggest to use not more than 48.
# Add export PATH=$PATH:/iplex/01/sys/applications/python/ibin to your .profile or .bashrc

#llrun -t 4 -c medium -g isimip -N Diffusion cpy runDiffusion.py
llrun -t 50 -c short -g isimip -N Diffusion cpy runDiffuse.py

