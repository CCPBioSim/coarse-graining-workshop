# CCPBioSim Coarse-Graining Workshop

[![ci](https://github.com/ccpbiosim/coarse-graining-workshop/actions/workflows/build.yaml/badge.svg?branch=main)](https://github.com/ccpbiosim/coarse-graining-workshop/actions/workflows/build.yaml)
[![latest](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fccpbiosim.github.io%2Fassets.json&query=%24.containers.coarse-graining-workshop.latest&labelColor=grey&logo=github&logoColor=white&label=latest&color=purple)](https://github.com/ccpbiosim/coarse-graining-workshop/pkgs/container/coarse-graining-workshop)
[![issues](https://img.shields.io/github/issues/ccpbiosim/coarse-graining-workshop?logo=github&labelColor=grey)](https://github.com/CCPBioSim/coarse-graining-workshop/issues)
[![pr](https://img.shields.io/github/issues-pr/ccpbiosim/coarse-graining-workshop?logo=github&labelColor=grey)](https://github.com/CCPBioSim/coarse-graining-workshop/pulls)

This workshop source repository contains the build recipe for a docker container derived from the CCPBioSim JupyterHub image. This container adds the necessary software packages and notebook content to form a deployable course container.

Here we present a practical introduction to the construction and simulation of coarse-grained (CG) membrane protein systems. We have also included a comparison to an equivalent all-atom (AA) simulation in order to raise awareness of where each has merit.

We have separated the content into 5 conceptually distinct parts:

1) Construction of an AA membrane protein system
2) Simulation/Visualisation of the AA system
3) Construction of a CG membrane protein system  
4) Equilibration of the CG system
5) Analysis of protein function and membrane dynamics

The user should work sequentially through each numbered directory within this repo\
starting from `1-allatom_setup` and finishing with `5-coarsegrain_analysis`


We acknowledge the work published by the `Vanni lab` which inspired the topic of this practical:

<b>Lipid scrambling is a general feature of protein insertases</b>\
<i>Li D., Rocha-Roa C., Schilling M.A., Vanni S.</i>\
https://www.pnas.org/doi/10.1073/pnas.2319476121

## How to Use

This training course is deployed on the [CCPBioSim](www.ccpbiosim.ac.uk) website via our cloud infrastructure, however you can deploy on your own machine with docker.

Pull the container from our repository::

    docker pull ghcr.io/ccpbiosim/coarse-graining-workshop:latest

In our containers we are using the JupyterHub default port 8888, so you should
forward this port when deploying locally::

    docker run -p 8888:8888 ghcr.io/ccpbiosim/coarse-graining-workshop:latest

## Authors

Workshop Content Authors:

- Robert Clark
- Iain Peter Shand Smith

## Contact

Please direct all questions and feedback to [Robert Clark](mailto:robert.clark@magd.ox.ac.uk) or [Iain Peter Shand Smith](mailto:iain.smith@bioch.ox.ac.uk)
