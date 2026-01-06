FROM python:3.10-slim

LABEL maintainer="BigBio Team <ypriverol@gmail.com>"
LABEL description="pyspectrafuse - Command-line utilities for spectrum clustering and conversion"

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY pyproject.toml .
COPY pyspectrafuse/ ./pyspectrafuse/

# Install the package
RUN pip install --no-cache-dir -e .

# Set the entrypoint
ENTRYPOINT ["pyspectrafuse_cli"]

# Default command
CMD ["--help"]

