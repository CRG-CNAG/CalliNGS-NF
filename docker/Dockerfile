FROM cbcrg/callings-nf@sha256:b65a7d721b9dd2da07d6bdd7f868b04039860f14fa514add975c59e68614c310

RUN apt-get update && apt-get install -y python \
 && curl -fsSL https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip > gatk-4.zip \
 && unzip gatk-4.zip \
 && rm gatk-4.zip

ENV PATH=/gatk-4.1.1.0:$PATH