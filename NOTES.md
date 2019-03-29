Files checked

NOTE: Scripts are used to develop final NF implementation (file)

Ok
-README.md - added GATK4 section at end
-vailidate-ci.sh - ok
-LICENSE - added NOTICE for GATK4 licence info
-/bin/gghist.R - ok
-/data - ok 
-/docker/Dockerfile - ok
-/figures - ok
-/scripts/4... - ok

Must Verify
-nextflow.config - line 6, calles GATK4 from DockerHub
-main.nf - main workflow, must run to verify breaking lines w/GATK4

Issues
-.travis.yml - line 15, pulls GATK docker image from AWS S3, fails requesting AWS credentials
-.circle.yml - line 10, pulls GATK docker image from AWS S3, fails requesting AWS credentials
-/docker/README.adoc - Line 15, points to GATK3 image (and Dockerfile - which includes GATK + other tools [Picard...])
-/scripts/1a... - line 19, comment points to GATK3
-/scripts/1b... - line 29, comment points to GATK3
-/scripts/2... - line 19 + 9 others, points to GATK3
-/scripts/3... - line 19 + 4 others, points to GATK3
-/scripts/5... - line 23 + 6 others, points to GATK3