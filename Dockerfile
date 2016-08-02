FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.



## Update Transform (should go away eventually)
#RUN \
#  . /kb/dev_container/user-env.sh && \
#  cd /kb/dev_container/modules && \
#  rm -rf transform && \ 
#  git clone https://github.com/kbase/transform -b develop
#
## setup the transform, but ignore errors because sample data cannot be found!
#RUN \
#  . /kb/dev_container/user-env.sh; \
#  cd /kb/dev_container/modules/transform/t/demo; \
#  python setup.py; \
#  exit 0;



# RUN apt-get update

# -----------------------------------------

# Install SDK Module
#
RUN mkdir -p /kb/module
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
WORKDIR /kb/module
RUN make


# Install MUSCLE
#
RUN mkdir -p /kb/module/muscle/bin
WORKDIR /kb/module/muscle/bin
RUN curl http://drive5.com/muscle/muscle3.8.425_binaries.tar.gz > muscle3.8.425_binaries.tar.gz
RUN tar xfz muscle3.8.425_binaries.tar.gz
RUN ln -s muscle3.8.425_i86linux64 muscle


WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
