FROM ubuntu:18.04

ENV DEPENDS='numpy==1.16 scipy==1.2 matplotlib==2.2' \
    dpr_version='0.1.2' \
    MPLBACKEND="agg"

RUN apt update && \
    apt install python3-pip -y --no-install-recommends && \
    apt autoclean && \
    # get python deps
    pip3 install setuptools wheel && \
    pip3 install $DEPENDS && \
    # install dpr itself
    pip3 install https://github.com/samuelstjean/dpr/releases/download/v${dpr_version}/dpr-${dpr_version}.tar.gz

# default command that will be run
CMD ["dpr","--help"]
