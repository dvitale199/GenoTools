FROM python:3.8.15

RUN apt update
RUN apt install git

COPY . /app

WORKDIR app

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN git clone https://github.com/dvitale199/GenoTools && cd GenoTools && pip install .

CMD ["python", "/app/run.py"]