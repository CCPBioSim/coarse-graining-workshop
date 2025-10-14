# Start with BioSim base image.
ARG BASE_IMAGE=latest
FROM ghcr.io/ccpbiosim/jupyterhub-base:$BASE_IMAGE

LABEL maintainer="James Gebbie-Rayet <james.gebbie@stfc.ac.uk>"

# Switch to jovyan user.
USER $NB_USER
WORKDIR $HOME

# Install workshop deps
RUN mamba install gromacs=2024.4=nompi_h5f56185_100 matplotlib=3.10.6 mdanalysis=2.9.0 numpy=2.3.3 scikit-learn=1.7.2 tqdm=4.67.1 mdtraj=1.11.0 iprogress ipywidgets salilab::dssp=3.0.0 -y
RUN pip install alphashape==1.3.1 click==8.3.0 click-log==0.4.0 coby==1.0.9 pbr==7.0.1 rtree==1.4.1 shapely==2.1.2 trimesh==4.8.3 vermouth==0.15.0

# Get workshop files and move them to jovyan directory.
COPY --chown=1000:100 . .
RUN rm -rf environment.yml LICENSE README.md docker .git .github

# Copy lab workspace
#COPY --chown=1000:100 docker/default-37a8.jupyterlab-workspace /home/jovyan/.jupyter/lab/workspaces/default-37a8.jupyterlab-workspace

# Always finish with non-root user as a precaution.
USER $NB_USER
