# ==========================
#  Ubuntu-based Pixi Image
# ==========================
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl ca-certificates \
    pkg-config libssl-dev libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="/root/.pixi/bin:$PATH"

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app

# Install environments with increased file descriptor limits
RUN bash -c "ulimit -n 4096 && pixi install -e default"
RUN bash -c "ulimit -n 4096 && pixi install -e analysis"

# Create the shell-hook bash script to activate the environment
RUN pixi shell-hook -e default -s bash > /shell-hook && \
    echo "#!/bin/bash" > /app/entrypoint.sh && \
    cat /shell-hook >> /app/entrypoint.sh && \
    echo 'exec "$@"' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

# Install fqkit in default environment
RUN pixi run install-fqkit

# Expose a port if your project serves something (optional)
EXPOSE 8000

# Entrypoint ensures environment is active before running any command
ENTRYPOINT ["/app/entrypoint.sh"]

# Default command (can be overridden when running container)
CMD ["bash"]
