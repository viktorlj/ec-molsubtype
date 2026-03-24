FROM python:3.12-slim

WORKDIR /app

# Install uv for fast dependency resolution
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Copy project files
COPY pyproject.toml README.md LICENSE ./
COPY src/ src/
COPY demo/ demo/

# Install the package
RUN uv pip install --system --no-cache .

# Cloud Run sets PORT env var (default 8080)
ENV PORT=8080
EXPOSE 8080

# Run with uvicorn, binding to 0.0.0.0:$PORT
CMD exec uvicorn ec_molsubtype.web.app:app --host 0.0.0.0 --port $PORT
