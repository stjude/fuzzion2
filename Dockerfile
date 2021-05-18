# syntax=docker/dockerfile:1

FROM debian:10 AS builder

ENV CPATH=/opt/htslib/include:$CPATH \
    LIBRARY_PATH=/opt/htslib/lib:$LIBRARY_PATH

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


FROM scratch as intermediate
WORKDIR /data

COPY --from=builder /lib/x86_64-linux-gnu/libc.so.6 /lib/x86_64-linux-gnu/libc.so.6
COPY --from=builder /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/libgcc_s.so.1
COPY --from=builder /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/libm.so.6
COPY --from=builder /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6
COPY --from=builder /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2

# ------ FUZZION2 -----
FROM intermediate as fuzzion2

COPY --from=builder /opt/htslib/lib/libhts.so.3 /lib/libhts.so.3
COPY --from=builder /lib/x86_64-linux-gnu/libpthread.so.0 /lib/x86_64-linux-gnu/libpthread.so.0
COPY --from=builder /usr/lib/x86_64-linux-gnu/libdeflate.so.0 /usr/lib/x86_64-linux-gnu/libdeflate.so.0
COPY --from=builder /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libz.so.1

COPY --from=builder /tmp/fuzzion2/build/bin/fuzzion2 /fuzzion2

ENTRYPOINT ["/fuzzion2"]
CMD ["-h"]

# ----- FUZZORT -----
FROM intermediate as fuzzort

COPY --from=builder /tmp/fuzzion2/build/bin/fuzzort /fuzzort

ENTRYPOINT ["/fuzzort"]

# ----- KMERANK -----
FROM intermediate as kmerank

COPY --from=builder /tmp/fuzzion2/build/bin/kmerank /kmerank

ENTRYPOINT ["/kmerank"]
CMD ["-h"]
