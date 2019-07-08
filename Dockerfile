FROM ubuntu:18.04

ENV DEPENDS='numpy==1.16.4 scipy==1.2.2 matplotlib==2.2.4' \
    dpr_version='0.1.2' \
    MPLBACKEND="agg"

RUN apt update && \
    apt install python3-pip -y --no-install-recommends && \
    apt autoclean && \
    # get python deps
    pip3 install --no-cache-dir setuptools wheel && \
    pip3 install --no-cache-dir $DEPENDS && \
    # install dpr itself
    pip3 install --no-cache-dir https://github.com/samuelstjean/dpr/releases/download/v${dpr_version}/dpr-${dpr_version}.tar.gz

# default command that will be run
CMD ["dpr","--help"]
