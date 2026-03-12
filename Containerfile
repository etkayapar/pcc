FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="5ecca81c4bb7f940b457e894c92bc96ec5b83cf9c6cf2ac1ef81977a7a8f97af"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/detect_outliers.yaml
#   prefix: /conda-envs/d802680873185ba15bf6619c7b7b3e9f
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - r-phytools =2.5_2
#       - r-ape =5.8_1
#       - coreutils
RUN mkdir -p /conda-envs/d802680873185ba15bf6619c7b7b3e9f
COPY envs/detect_outliers.yaml /conda-envs/d802680873185ba15bf6619c7b7b3e9f/environment.yaml

# Conda environment:
#   source: envs/iqtree.yaml
#   prefix: /conda-envs/a2e8cf3555187760f541735df9480b1e
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - iqtree =2.3
#       - coreutils
RUN mkdir -p /conda-envs/a2e8cf3555187760f541735df9480b1e
COPY envs/iqtree.yaml /conda-envs/a2e8cf3555187760f541735df9480b1e/environment.yaml

# Conda environment:
#   source: envs/mafft.yaml
#   prefix: /conda-envs/6c2b0c8fa0679fb27d1081bae6230c9e
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - mafft =7.525
#       - trimal =1.5.1
#       - biopython =1.83
#       - coreutils
RUN mkdir -p /conda-envs/6c2b0c8fa0679fb27d1081bae6230c9e
COPY envs/mafft.yaml /conda-envs/6c2b0c8fa0679fb27d1081bae6230c9e/environment.yaml

# Conda environment:
#   source: envs/phylo_scripts_python.yaml
#   prefix: /conda-envs/0db6a624dc70e8008365f5f11a06893a
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - biopython =1.83
#       - coreutils
RUN mkdir -p /conda-envs/0db6a624dc70e8008365f5f11a06893a
COPY envs/phylo_scripts_python.yaml /conda-envs/0db6a624dc70e8008365f5f11a06893a/environment.yaml

# Conda environment:
#   source: envs/treeshrink.yaml
#   prefix: /conda-envs/ab9381e9b1804d68d3cc2ab42f167fc6
#   name: treeshrink
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - treeshrink =1.3.9
#     - coreutils
RUN mkdir -p /conda-envs/ab9381e9b1804d68d3cc2ab42f167fc6
COPY envs/treeshrink.yaml /conda-envs/ab9381e9b1804d68d3cc2ab42f167fc6/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/d802680873185ba15bf6619c7b7b3e9f --file /conda-envs/d802680873185ba15bf6619c7b7b3e9f/environment.yaml && \
    conda env create --prefix /conda-envs/a2e8cf3555187760f541735df9480b1e --file /conda-envs/a2e8cf3555187760f541735df9480b1e/environment.yaml && \
    conda env create --prefix /conda-envs/6c2b0c8fa0679fb27d1081bae6230c9e --file /conda-envs/6c2b0c8fa0679fb27d1081bae6230c9e/environment.yaml && \
    conda env create --prefix /conda-envs/0db6a624dc70e8008365f5f11a06893a --file /conda-envs/0db6a624dc70e8008365f5f11a06893a/environment.yaml && \
    conda env create --prefix /conda-envs/ab9381e9b1804d68d3cc2ab42f167fc6 --file /conda-envs/ab9381e9b1804d68d3cc2ab42f167fc6/environment.yaml && \
    conda clean --all -y
