# oesphlora
# oesphlora
![CI](https://github.com/mauricebarrett/oesphlora/actions/workflows/ci.yml/badge.svg?branch=main)
![Release](https://github.com/mauricebarrett/oesphlora/actions/workflows/release.yml/badge.svg?branch=main)
[![semantic-release](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--release-e10079.svg)](https://github.com/semantic-release/semantic-release)
![Latest Release](https://img.shields.io/github/v/release/mauricebarrett/oesphlora?sort=semver)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-Angular-brightgreen.svg)](https://www.conventionalcommits.org/en/v1.0.0/)


## Overview
This repository contains code and analysis scripts accompanying the manuscript entitled:
"Alterations in the oesophago-gastric mucosal microbiome in patients along the inflammation-metaplasia-dysplasia-oesophageal adenocarcinoma sequence"
Currently under preparation.


## Installation

There are two ways to set up this project:

1. Docker
2. Native installation (Only for Ubuntu)

### Option 1: Docker

The easiest way is to use the prebuilt Docker image, which will contain required dependencies.

### Option 2: Native installation

This project uses [pixi](https://pixi.sh) for environment management.
Most dependencies are handled automatically by pixi.

#### Step 1: Install pixi maneged dependacies
Follow the instuctions in the link to install pixi [pixi installation guide](https://pixi.sh/latest/#installation)

Once pixi is install one may enter the

```bash
pixi intall
```

This should install all dependacies installed by

#### Step 2: Install all other dependacies

⚠️ **Important:** Some tools are a system dependency and **cannot be installed through pixi**.
You must install it separately before using reproducing analysis


##### Install `fqkit`



## Author

**Maurice Barrett**  
University College Cork  
Email: mauricepatrickbarrett@gmail.com  
