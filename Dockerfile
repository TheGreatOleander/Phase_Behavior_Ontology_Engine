FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install dependencies first (cached layer)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY *.py ./

# Copy frontend (index.html lives at project root; templates/ and static/ are
# optional — we use a conditional so the build succeeds even if absent)
COPY index.html ./
RUN mkdir -p /app/templates /app/static

# Create data directory for optional topological DB
RUN mkdir -p /app/data

# Non-root user for security
RUN useradd -m appuser && chown -R appuser /app
USER appuser

# Environment
ENV FLASK_APP=app.py
ENV FLASK_DEBUG=false
ENV PORT=5000

EXPOSE 5000

# Use gunicorn for production (not Flask's dev server)
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "2", "--timeout", "120", "app:app"]
