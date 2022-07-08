FROM python:3.8-slim-bullseye

RUN apt update && apt install -y git

RUN git clone https://github.com/cov-lineages/constellations.git && \
    git clone https://github.com/GenomePathogenAnalysisService/gpas-testing.git

WORKDIR /gpas-testing
RUN git clone https://github.com/oxfordmmm/gumpy.git

WORKDIR /gpas-testing
RUN pip3 install .
