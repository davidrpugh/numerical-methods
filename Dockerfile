FROM andrewosh/binder-base

MAINTAINER davidrpugh <david.pugh@maths.ox.ac.uk>

USER root

RUN apt-get update && \
    apt-get install -y  --no-install-recommends dynare && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER main
