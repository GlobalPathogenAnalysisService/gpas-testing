FROM python:3.8-slim-bullseye

RUN apt update && apt install -y git

RUN git clone https://github.com/phe-genomics/variant_definitions.git && \
    git clone https://github.com/cov-lineages/constellations.git && \
    git clone https://github.com/GenomePathogenAnalysisService/gpas-covid-synthetic-reads.git

WORKDIR /gpas-covid-synthetic-reads
RUN git clone https://github.com/oxfordmmm/gumpy.git

WORKDIR /gpas-covid-synthetic-reads/gumpy
RUN python3 -m pip install --upgrade pip && \
    pip3 install -r requirements.txt && \
    python3 setup.py install --user

WORKDIR /gpas-covid-synthetic-reads
RUN pip3 install -r requirements.txt && \
    python3 setup.py install --user
