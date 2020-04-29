# sudo docker pull debian:buster
FROM debian:buster

RUN apt-get update && apt-get install -y \
    ncbi-entrez-direct \
    ncbi-blast+ \
    aragorn \
    python3-pip \
    libdatetime-perl \
    libxml-simple-perl \
    libdigest-md5-perl \
    wget \
    default-jre-headless \
    bioperl

RUN pip3 install click biopython bcbio-gff prefect

# install prokka
#RUN cpan Bio::Perl
RUN wget -O /prokka.zip https://github.com/tseemann/prokka/archive/v1.14.6.zip; \
    unzip -d / /prokka.zip; rm /prokka.zip; \
    # remove some unused databases
    rm -r /prokka-1.14.6/db/kingdom/Archaea \
    /prokka-1.14.6/db/kingdom/Mitochondria \
    /prokka-1.14.6/db/kingdom/Viruses \
    /prokka-1.14.6/db/genus/*; \
    # consider including Pfam-A when the new release comes
    /prokka-1.14.6/bin/prokka --setupdb; \
    rm -r /prokka-1.14.6/db/cm/Archaea* /prokka-1.14.6/db/cm/Viruses*
ENV PATH="$PATH:/prokka-1.14.6/bin"

# install PipolinFinder
COPY ExplorePipolin-master.zip .
# replace COPY with wget when the repo is public
RUN unzip -d / /ExplorePipolin-master.zip; rm /ExplorePipolin-master.zip; \
    pip3 install /ExplorePipolin-master; rm -r /ExplorePipolin-master

ENTRYPOINT ["pipolin_finder"]

