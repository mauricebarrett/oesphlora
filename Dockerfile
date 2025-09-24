# ==========================
#  Stage 1: Build with Pixi
# ==========================
FROM ghcr.io/prefix-dev/pixi:0.40.0 AS build

# Copy project files into /app
COPY . /app
WORKDIR /app

# Install all environments
RUN pixi install -e default
RUN pixi install -e setup
RUN pixi install -e data-processing


# Run setup task (installs fqkit into /app/.cargo/bin)
RUN pixi run -e setup setup


# Generate a shell hook that sets up the environment
RUN pixi shell-hook -e default > /shell-hook.sh

# Add a line so any command passed to "docker run" executes inside the env
RUN echo 'exec "$@"' >> /shell-hook.sh


# ==========================
#  Stage 2: Production image
# ==========================
FROM ubuntu:24.04 AS production

# Set working directory
WORKDIR /app

# Copy environment + project code + shell hook from build stage
COPY --from=build /app /app
COPY --from=build /shell-hook.sh /shell-hook.sh

# Expose a port if your project serves something (optional)
EXPOSE 8000

# Entrypoint ensures environment is active before running any command
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]

# Default command (can be overridden when running container)
CMD ["bash"]
