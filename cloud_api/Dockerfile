# Base image with essential dependencies
FROM python:3.10.13

# Set working directory
WORKDIR /app

# Copy requirements and project files
COPY requirements.txt /app
COPY . /app

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Install the_real_genotools specifically
RUN pip install the_real_genotools

# Expose port for Streamlit
EXPOSE 8501

# Start Streamlit app
CMD ["streamlit", "run", "app.py"]