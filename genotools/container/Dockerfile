FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    gcc \
    build-essential \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY . /app

WORKDIR app

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN git clone https://github.com/dvitale199/GenoTools && cd GenoTools && pip install .

CMD ["python", "/app/run.py"]