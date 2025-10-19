FROM python:3.12-slim as builder

# Install node (for building the frontend)
RUN apt-get update && apt-get install -y curl gnupg ca-certificates build-essential && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash -
RUN apt-get update && apt-get install -y nodejs && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy frontend and build it
COPY frontend/package.json frontend/package-lock.json* frontend/
RUN cd frontend && npm install && npm run build

# Install Python deps
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app source
COPY . .

FROM python:3.12-slim
WORKDIR /app

# Copy built frontend and installed python packages from builder
COPY --from=builder /app /app

ENV PYTHONUNBUFFERED=1

EXPOSE 8000

CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
