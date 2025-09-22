# ==========================
#  Stage 1: Build with Pixi
# ==========================
FROM ghcr.io/prefix-dev/pixi:0.40.0 AS build

# Copy project files into /app
COPY . /app
WORKDIR /app

# Install dependencies into the "prod" environment
# (make sure you have an [env.prod] defined in pixi.toml)
RUN pixi install -e prod

# Generate a shell hook that sets up the environment
RUN pixi shell-hook -e prod > /shell-hook.sh

# Add a line so any command passed to "docker run" executes inside the env
RUN echo 'exec "$@"' >> /shell-hook.sh


# ==========================
#  Stage 2: Production image
# ==========================
FROM ubuntu:24.04 AS production

# Set working directory
WORKDIR /app

# Copy environment + project code + shell hook from build stage
COPY --from=build /app/.pixi/envs/prod /app/.pixi/envs/prod
COPY --from=build /app /app
COPY --from=build /shell-hook.sh /shell-hook.sh

# Expose a port if your project serves something (optional)
EXPOSE 8000

# Entrypoint ensures environment is active before running any command
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]

# Default command (can be overridden when running container)
CMD ["bash"]
