FROM debian:buster-slim

RUN apt-get update && \
  apt-get install -y --no-install-recommends \
    autoconf \
    autoconf-archive \
    automake \
    build-essential \
    git \
    libcppunit-dev \
    valgrind && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/ && \
  rm -rf /tmp/downloaded_packages/ /tmp/*.rds 

RUN adduser developer

USER developer
