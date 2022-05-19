FROM python:3.8-slim-bullseye

RUN apt update && apt install -y git

RUN git clone https://github.com/cov-lineages/constellations.git && \
    git clone https://github.com/GenomePathogenAnalysisService/gpas-covid-synthetic-reads.git

WORKDIR /gpas-covid-synthetic-reads
RUN git clone https://github.com/oxfordmmm/gumpy.git

WORKDIR /gpas-covid-synthetic-reads/gumpy
RUN python3 -m pip install --upgrade pip && \
    pip3 install -r requirements.txt && \
    python setup.py build --force && \
    pip3 install .

WORKDIR /gpas-covid-synthetic-reads
RUN pip3 install -r requirements.txt && \
    pip3 install .
