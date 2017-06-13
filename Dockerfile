FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.


# RUN apt-get update

# -----------------------------------------
# Install MUSCLE
#
RUN mkdir -p /kb/module/muscle/bin && \
    cd /kb/module/muscle/bin && \
    curl http://drive5.com/muscle/muscle3.8.425_binaries.tar.gz > muscle3.8.425_binaries.tar.gz && \
    tar xfz muscle3.8.425_binaries.tar.gz && \
    ln -s muscle3.8.425_i86linux64 muscle && \
    rm muscle3.8.425_binaries.tar.gz 

# Install SDK Module
#
#RUN mkdir -p /kb/module
ADD ./ /kb/module
WORKDIR /kb/module
RUN \
    mkdir -p /kb/module/work && \
    make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
