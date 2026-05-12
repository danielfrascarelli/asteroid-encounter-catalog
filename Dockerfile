FROM python:3.11-slim

# Build deps for compiled extensions (scipy, rebound)
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Install uv from official image (pinned major version)
COPY --from=ghcr.io/astral-sh/uv:0.4 /uv /usr/local/bin/uv

# Install directly to system Python — no .venv management inside container
ENV UV_SYSTEM_PYTHON=1

WORKDIR /app

# Dependency layer (cached unless pyproject.toml changes)
COPY pyproject.toml ./
RUN uv pip install -e ".[dev]" --no-cache

# Source code
COPY . .

# data/ and logs/ are always mounted from the host (see docker-compose.yml)
VOLUME ["/app/data", "/app/logs"]

CMD ["python", "-m", "scripts.run_pipeline", "--config", "config.yaml"]
