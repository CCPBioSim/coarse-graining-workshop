# Start with BioSim base image.
ARG BASE_IMAGE=latest@sha256:5eb67b09e3fdf64b3128bd9501e08a0cbb007804b678f2096728d8c1260e6217
FROM ghcr.io/ccpbiosim/jupyterhub-base:$BASE_IMAGE

LABEL maintainer="James Gebbie-Rayet <james.gebbie@stfc.ac.uk>"

ARG TARGETPLATFORM

# Switch to jovyan user.
USER $NB_USER
WORKDIR $HOME

# Install workshop deps
RUN conda install -y -c conda-forge \
    conda-forge::matplotlib=3.11.2 \
    conda-forge::mdanalysis=2.10.0 \
    conda-forge::numpy=2.5.3 \
    conda-forge::scikit-learn=1.9.1 \
    conda-forge::tqdm=4.67.1 \
    conda-forge::mdtraj=1.11.1 \
    conda-forge::iprogress=0.4 \
    conda-forge::ipywidgets=8.1.9 \
    conda-forge::gromacs=2026.3=nompi_*
RUN pip install --no-cache-dir \
    alphashape==1.3.1 \
    click==8.5.0 \
    click-log==0.4.0 \
    coby==1.0.15 \
    pbr==7.0.3 \
    rtree==1.4.1 \
    shapely==2.1.2 \
    trimesh==4.12.2 \
    vermouth==0.15.0
    
# Get workshop files and move them to jovyan directory.
COPY --chown=1000:100 . .
RUN rm -rf environment.yml LICENSE README.md docker .git .github

# Copy lab workspace
#COPY --chown=1000:100 docker/default-37a8.jupyterlab-workspace /home/jovyan/.jupyter/lab/workspaces/default-37a8.jupyterlab-workspace

# Always finish with non-root user as a precaution.
USER $NB_USER
