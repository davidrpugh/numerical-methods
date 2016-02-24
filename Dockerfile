FROM andrewosh/binder-base

MAINTAINER davidrpugh <david.pugh@maths.ox.ac.uk>

USER root

# Install Dynare
RUN apt-get update && \
    apt-get install -y  --no-install-recommends dynare && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install additional Python dependencies
RUN conda update conda && \
    conda install -y seaborn

USER main
