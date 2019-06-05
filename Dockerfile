FROM kbase/sdkbase2:python
MAINTAINER KBase Developer [Dylan Chivian (DCChivian@lbl.gov)]
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Install MUSCLE
#
RUN mkdir -p /kb/module/muscle/bin && \
    cd /kb/module/muscle/bin && \
    curl http://drive5.com/muscle/muscle3.8.425_binaries.tar.gz > muscle3.8.425_binaries.tar.gz && \
    tar xfz muscle3.8.425_binaries.tar.gz && \
    ln -s muscle3.8.425_i86linux64 muscle && \
    rm muscle3.8.425_binaries.tar.gz 

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
