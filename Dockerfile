# syntax=docker/dockerfile:1

FROM debian:10 AS builder

ENV CPATH=/opt/htslib/include:$CPATH \
    LIBRARY_PATH=/opt/htslib/lib:$LIBRARY_PATH

# fuzzion2: g++ make
# htslib: gcc lbzip2 libdeflate-dev make wget zlib1g-dev
# wget: ca-certificates
RUN apt-get update \
        && apt-get --yes install --no-install-recommends \
            ca-certificates \
            g++ \
            gcc \
            lbzip2 \
            libdeflate-dev \
            make \
            wget \
            zlib1g-dev \
        && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# This build of htslib disables the bzip2 and LZMA CRAM block compression
# methods and network support.
RUN wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 \
        && echo "2280141b46e953ba4ae01b98335a84f8e6ccbdb6d5cdbab7f70ee4f7e3b6f4ca  htslib-1.12.tar.bz2" | sha256sum --check \
        && tar xf htslib-1.12.tar.bz2 \
        && cd htslib-1.12 \
        && ./configure --prefix /opt/htslib --disable-bz2 --disable-libcurl --disable-lzma \
        && make --jobs $(ncpus) \
        && make install \
        && rm -r /tmp/htslib-1.12*

WORKDIR /tmp/fuzzion2

COPY Makefile ./
COPY src/ src/

RUN make --jobs $(ncpus)

FROM debian:10

ENV PATH=/opt/fuzzion2/bin:$PATH \
    LD_LIBRARY_PATH=/opt/fuzzion2/lib:/opt/htslib/lib:$LD_LIBRARY_PATH

COPY --from=builder /opt/htslib/lib/ /opt/htslib/lib/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libdeflate.so.0 /opt/fuzzion2/lib/
COPY --from=builder /tmp/fuzzion2/build/bin/ /opt/fuzzion2/bin/

CMD ["fuzzion2"]
