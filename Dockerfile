FROM python:2
LABEL maintainer="Uemit Seren <uemit.seren@gmail.com>"

RUN useradd --uid=10372 -ms /bin/bash gwas-web

RUN mkdir /srv/data
VOLUME '/srv/data'

WORKDIR /gwapp
COPY backend/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN chown -R gwas-web:root /gwapp && chmod -R 755 /gwapp

RUN echo "cm_dir = /srv/data/gwas-web/dataset/" > /home/gwas-web/.gwa_config

# required because of https://github.com/docker/docker/issues/20240
RUN cd backend && ln -s static ../frontend/target/gwaswebapp-0.0.1-SNAPSHOT/ && cd public && ln -s resources ../../frontend/src/main/java/com/gmi/gwaswebapp/client/resources/ && cd /gwapp

ENV PYTHONPATH=/gwapp/atgwas/src/:/gwapp/:/gwapp/atpipeline/
ENV GOTO_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
USER gwas-web

EXPOSE 8080

ENTRYPOINT ["python" , "backend/__init__.py"]