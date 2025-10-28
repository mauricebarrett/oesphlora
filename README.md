# oesphlora

![CI](https://github.com/mauricebarrett/oesphlora/actions/workflows/ci.yml/badge.svg?branch=main)
![Release](https://github.com/mauricebarrett/oesphlora/actions/workflows/release.yml/badge.svg?branch=main)
[![semantic-release: angular](https://img.shields.io/badge/semantic--release-angular-e10079.svg?logo=semantic-release)](https://github.com/semantic-release/semantic-release)
[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)
[![GitHub release](https://img.shields.io/github/v/release/OWNER/REPO)](https://github.com/OWNER/REPO/releases)
[![Docker Pulls](https://img.shields.io/docker/pulls/DOCKERHUB_USER/IMAGE)](https://hub.docker.com/r/DOCKERHUB_USER/IMAGE)


## Overview
This repository contains code and analysis scripts accompanying the manuscript entitled:
"Alterations in the oesophago-gastric mucosal microbiome in patients along the inflammation-metaplasia-dysplasia-oesophageal adenocarcinoma sequence"
Currently under preparation.


## Installation

There are two ways to set up this project:

1. Docker
2. Native installation (Only for Ubuntu)

either way you will need to clone the repo

```bash
git clone https://github.com/mauricebarrett/oesphlora.git
cd oesphlora
```

### Option 1: Docker

The easiest way is to use the prebuilt Docker image, which will contain required dependencies.

#### Step 1: Install Docker

Follow the instuctions in the link to install Docker [Install Docker Engine] https://docs.docker.com/engine/install/

### Option 2: Native installation

This project uses [pixi](https://pixi.sh) for environment management.
Most dependencies are handled automatically by pixi.

#### Step 1: Install Pixi
Follow the instuctions in the link to install pixi [pixi installation guide](https://pixi.sh/latest/#installation)

#### Step 2: Install the deafult enviroment
Install deafult enviroment as follows

```bash
pixi intall -e default
```

#### Step 3: Installing dependacies not manged by Conda (fqkit)
Some dependencies cannot be installed via Pixi because they are not managed by Conda. This workflow uses fqkit, which is installed from its Rust crate using cargo. The command below installs fqkit into a local cargo directory within the project.

```bash
pixi run install-fqkit
```

#### Step 4: Install main analysis enviroment
Install main analysis enviroment


## Author

**Maurice Barrett**
University College Cork
Email: mauricepatrickbarrett@gmail.com
