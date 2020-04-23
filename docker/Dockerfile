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
RUN wget -O $HOME/prokka.zip https://github.com/tseemann/prokka/archive/v1.14.6.zip; \
unzip -d $HOME $HOME/prokka.zip; rm $HOME/prokka.zip; \
# remove unused databases
rm -r $HOME/prokka-1.14.6/db/kingdom/[Archaea,Mitochondria,Viruses] \
$HOME/prokka-1.14.6/db/genus/* $HOME/prokka-1.14.6/db/cm/[Archaea,Viruses]
# consider including Pfam-A when the new release come
RUN $HOME/prokka-1.14.6/bin/prokka --setupdb

