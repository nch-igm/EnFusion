FROM rocker/r-ver:3.5.0
#Use rocker from well maininted dockerhub as the parent image, will set up as a debian system
#debian:stretch

LABEL software.name="EnFusion"
LABEL software.version="EnFusion v1.0.0"
LABEL software.description="ensemble approach for fusion detection from RNA-seq data"
LABEL container.base.image="debian"
LABEL tags="Fusion Calling"

#Install necessary tools onto debian system
run apt-get update
run apt-get --yes install vim
run apt-get --yes install less
run apt-get --yes install parallel
run apt-get --yes install python-pip
run apt-get --yes install curl
run apt-get --yes install libpq-dev

#Install necessary R packages used in overlap script
RUN R -e "install.packages('devtools', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('optparse', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('VennDiagram', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('gridExtra', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('dplyr', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('readr', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('DBI', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('RPostgres', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('dbplyr', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('magrittr', repos='http://cran.us.r-project.org/')"
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org/')"

#Install tools using pip that will be used in the analysis
run pip install awscli
run pip install boto3
run apt-get --yes install jq


COPY SCRIPTS /SCRIPTS

#ARG s3_folder=

# Run
ENTRYPOINT ["bin/bash",  "/SCRIPTS/kickoff_overlap.sh"]

#CMD [ "--output_location", "test_data" ]
