FROM ghcr.io/scientificcomputing/fenics-gmsh:2023-11-15

# We set the working directory to install docker dependencies
WORKDIR /tmp/

# Install the Python-module
COPY . /tmp
RUN python3 -m pip install  --no-cache-dir .


# We remove the contents of the temporary directory to minimize the size of the image
RUN rm -rf /tmp
