FROM python:3.11-bookworm

RUN apt-get update
# singularity build prerequisites
RUN apt-get install -y \
   autoconf \
   automake \
   cryptsetup \
   git \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev \
   make \
   cmake \
   sudo \
   virtualenv \
   tree

WORKDIR /root

RUN virtualenv venv
ENV py3=/root/venv/bin/python3
RUN $py3 -m pip install --upgrade pip
RUN $py3 -m pip install spython


# install go, which singularity is written in
COPY installGO.sh installGO.sh
RUN sh installGO.sh

ENV PATH="/usr/local/go/bin:$PATH"

COPY installSingularity.sh installSingularity.sh
ENV SINGULARITY_VERSION=4.1.0
RUN sh installSingularity.sh
	
ENV PATH="/root/local/singularity-ce-${SINGULARITY_VERSION}:$PATH"


COPY runsingularity.sh /root/runsingularity.sh
RUN chmod +x /root/runsingularity.sh

# this needs to be run with elevated privileges (hence the --privileged flag in the docker run command)
ENTRYPOINT ["/root/runsingularity.sh"]
